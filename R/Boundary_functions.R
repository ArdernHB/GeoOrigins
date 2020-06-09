


#' Returns geographic boundaries identified from trait distances
#'
#' This function takes the distances among all reference specimens and uses a process similar to
#' leave-one-out correct cross-validation to identify likely spatial trait boundaries. This method
#' only allows for the inclusion of specimens with specific known locations (i.e. specific
#' latitude-longitude coordinates). This functions requires a distance input, which allows the user
#' to input any desired dissimilarity or distance metrics as appropriate for the data.
#' @param RefDistMat is a square matrix of pairwise distances among all reference specimens.
#' @param ExpandMap is a vector of 2 elements for expanding the plotting region of the map. The first element expands the latitudinal area and the second element expands the longitudinal area.
#' @param DataDump determines whether a datadump is used, which is recommended for large datasets that may take time to computeand when the user might want to interupt the process mid way. The default is set to TRUE. If it's set to FALSE then the results are output as an array where rows and columns reflect latitude longitude and the third dimension are individual specimens in the reference dataset.
#' @param DataDumpPath is a directory path indicating where the data for constructing the boundaries can be dumped.
#' @param StartPoint is an iteger that denotes the specimen to start the process on. As the method needs to cycle through all known specimens in the database, this can take some time. If the process is stopped for whatever reason the process can be picked up again by adjusting the \code{StartPoint} to the specimen number that it had previously finished on.
#' @param RefIDs is a vector of the unique identifiers for each of the reference specimens in the reference dataset. These values will be used for naming the files in the datadump folder and also for matching up with the returned summary results. Default is set to NULL and if this is not populated then reference data is worked through consecutively, naming the datadump files in consecutive order.
#' @param IgnorePrompts default is set to FALSE, but if set to TRUE queries such as those confirming the location of the datadump will be suppressed.
#' @inheritParams IDbyDistanceDistInput
#' @return If Verbose is TRUE then an array of all values for every sampled grid reference is returned. If Verbose is FALSE then only those grid references with the highest correlation values are returned.
#' @details The map plotting of this function makes use of the functions of the \code{maps} package.
#' Spatial polygons were constructed using functions from the \code{sp} package and calculations of polygon intersects and areas were carried out using functions from the \code{raster} package.
#' When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
#'
#'
#' @section Citations:
#'
#' Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
#' Enhancements by Thomas P Minka and Alex Deckmyn. (2017). maps: Draw Geographical Maps. R
#' package version 3.2.0. https://CRAN.R-project.org/package=maps
#'
#' Grosjean, Ph. (2016). svMisc: SciViews-R. UMONS, Mons, Belgium.
#' http://www.sciviews.org/SciViews-R.
#'
#' Pebesma, E.J., R.S. Bivand, 2005. sp: Classes and methods for spatial data in R. R News 5 (2)
#' http://cran.r-project.org/doc/Rnews/.
#'
#' Robert J. Hijmans (2016). raster: Geographic Data Analysis and Modeling. R package version 2.5-8.
#' https://CRAN.R-project.org/package=raster
#'
#' @author Ardern Hulme-Beaman
#'
#' @examples
#' RatDistMat <- ProcDistanceTable(Rpraetor$LMs)
#'
#' Range.Exp <- .5
#'
#' Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Long)+Range.Exp))
#' Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Lat)+Range.Exp))
#'
#' rThres <- IDbyDistanceDistInputCCV(LatLongs = Rpraetor$Lat.Long,
#'               DistDataMat = RatDistMat,
#'               Verbose = TRUE,
#'               ProvConfidence = .95,
#'               PrintProg = FALSE,
#'               Method = 'Spearman')
#'
#' LowResRsamp.Rp <- c(7,20)
#' Boundaryfinding <- BoundaryFinder(LatLongs = Rpraetor$Lat.Long,
#'                                RefDistMat = RatDistMat,
#'                                LongRange = Long.Range,
#'                                LatRange = Lat.Range,
#'                                RangeSamp = LowResRsamp.Rp,
#'                                PlotValCor = rThres$`Provenancing.Correlation.95%.Confidence`,
#'                                ExpandMap = c(0,0),
#'                                RefIDs = rownames(Rpraetor$Lat.Long),
#'                                DataDump = FALSE,
#'                                IgnorePrompts = TRUE,
#'                                Method = 'Spearman')
#'
#' @import rgeos
#' @import maps
#' @import svMisc
#' @export



BoundaryFinder <- function(LatLongs, RefDistMat=matrix(), LongRange, LatRange, RangeSamp=10, PrintProg=TRUE, PlotValCor, ExpandMap=c(0,0), DataDump=TRUE, DataDumpPath=NA, StartPoint=1, RefIDs=NULL, IgnorePrompts=FALSE, Method=c('Pearson', 'Spearman')){

  UserInputAssessment(LatLongs, RefDistMat, LongRange, LatRange, RangeSamp, Method, RefData = NA, DistVec = NA)


  if (IgnorePrompts==FALSE){
    if (DataDump==TRUE){
      if (is.na(DataDumpPath)){
        stop("\n Datadump directory not set. Please set DataDumpPath argument and run function again.")
      } else {
        if ((length(list.files(DataDumpPath))+length(list.dirs(DataDumpPath)))>1){
          QueryDD <- paste("\n Datadump directory",  DataDumpPath, " contains existing files.\n ",
                           "Do you wish to continue with this directory? \n please respond `y` or `n` \n \n", sep = '')

          UserResponse <- readline(QueryDD)


          UserGood <- 'Bad'
          while (UserGood!="Good"){
            if (!(UserResponse %in% c('y','n'))){
              UserResponse <- readline("User response was neither y or n. Please respond either y or n (you are stuck in a while loop!)")
            } else if (UserResponse == 'n'){
              stop("User response is n. Please change the DataDumpPath argument to required location and rerun function")
            } else if (UserResponse == 'y'){
              print(paste("User response was y. Function will proceed and save files to ", DataDumpPath, sep=''))
              UserGood <- "Good"
            }
          }


        } else {
          QueryDD <- paste("\n DataDumpPath directory set to",  DataDumpPath,
                           "Do you wish to continue with this directory? \n please respond y or n \n \n",
                           sep = '')

          UserResponse <- readline(QueryDD)


          UserGood <- 'Bad'
          while (UserGood!="Good"){
            if (!(UserResponse %in% c('y','n'))){
              UserResponse <- readline("User response was neither y or n. Please respond either y or n (you are stuck in a while loop!)")
            } else if (UserResponse == 'n'){
              stop("User response is n. Please change the DataDumpPath argument to required location and rerun function")
            } else if (UserResponse == 'y'){
              print(paste("User response was y. Function will proceed and save files to ", DataDumpPath, sep=''))
              UserGood <- "Good"
            }
          }


        }

      }
    } else {
      QueryDD <- paste("\n Datadump is set to FALSE, so data will be stored as an array and returned at the end of the process.\n If process is halted all data will be lost. \n ",
                       "Do you wish to continue with this and have data stored locally? \n please response `y` or `n` \n \n",
                        sep = '')

      UserResponse <- readline(QueryDD)


      UserGood <- 'Bad'
      while (UserGood!="Good"){
        if (!(UserResponse %in% c('y','n'))){
          UserResponse <- readline("User response was neither y or n (you are stuck in a while loop! you must respond correctly). Please respond either y or n")
        } else if (UserResponse == 'n'){
          stop("User response is n. Please change the DataDump argument and set required location via DataDumpPath argument and rerun function")
        } else if (UserResponse == 'y'){
          print(paste("User response was y. Function will proceed and store results as an array to be returned on completion of process", sep=''))
          warning("\n WARNING DataDump was set to FALSE. Results are saved to an array and if the process was stopped all data will have been lost.")
          UserGood <- "Good"
        }
      }

    }
  }



  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")

  #creating a range that will cover the whole geographic area of interest
  #this is for the function to loop through later on
  if (length(RangeSamp)==1){
    LongSamp <- RangeSamp
    LatSamp <- RangeSamp
  } else if (length(RangeSamp)==2){
    LongSamp <- RangeSamp[2]
    LatSamp <- RangeSamp[1]
  } else if (length(RangeSamp)>2){
    stop("too many dimensions in RangeSamp")
  }

  LongRangeSteps <- (LongRange[2]-LongRange[1])/(LongSamp-2)
  LatRangeSteps <- (LatRange[2]-LatRange[1])/(LatSamp-2)

  #this output for Lat/Long ways provides what the loop should sequence through
  Longways <- c(LongRange[1]-LongRangeSteps, seq(LongRange[1], LongRange[2], by = LongRangeSteps), LongRange[2]+LongRangeSteps)
  Latways <- c(LatRange[1]-LatRangeSteps, seq(LatRange[1], LatRange[2], by = LatRangeSteps), LatRange[2]+LatRangeSteps)

  #plotting a map that will later be populated with polygons for each jackknifing iteration in the dd loop below
  if (sum(ExpandMap>0)>0){
    maps::map("world", xlim=c(min(Longways)-ExpandMap[2], max(Longways)+ExpandMap[2]), ylim=c(min(Latways)-ExpandMap[1], max(Latways)+ExpandMap[1]), interior=FALSE, col="black", bg=graphics::par(bg="white"))
    graphics::points(x = as.character(LatLongs$Long), y = as.character(LatLongs$Lat), pch=23, bg='blue',  cex=1)
  } else {
    maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"))
    graphics::points(x = as.character(LatLongs$Long), y = as.character(LatLongs$Lat), pch=23, bg='blue',  cex=1)
  }

  SpecimenLoc <- cbind(chr2nu(LatLongs$Long), chr2nu(LatLongs$Lat))
  DistPol <-  sp::Polygon(SpecimenLoc[grDevices::chull(SpecimenLoc),])
  DistPols <-  sp::Polygons(list(DistPol),1)
  DistSpatPol <-  sp::SpatialPolygons(list(DistPols))
  sp::proj4string(DistSpatPol) <- "+proj=longlat +ellps=WGS84"
  DistSpatPol$area <- raster::area(DistSpatPol)


  if (is.null(RefIDs)){
    RefIDs <- 1:dim(RefDistMat)[1]
  }

  if (sum(rownames(RefDistMat)==c(1:length(rownames(RefDistMat))))==length(rownames(RefDistMat))){
    ArrayDimNames <- list(Latways, Longways, rownames(RefDistMat))
    Prov.Array.Res <- array(NA, dim = c(length(Latways), length(Longways), length(StartPoint:dim(RefDistMat)[1])), dimnames = ArrayDimNames)

    if (IgnorePrompts==FALSE){
      print(x = "Distance matrix rownames appear match numeric order, so order number has been used as IDs for DataDump file names or array dimensions.")
    }
  } else if (length(unique(rownames(RefDistMat)))==length(rownames(RefDistMat))){
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(RefDistMat)))
    Prov.Array.Res <- array(NA, dim = c(length(Latways), length(Longways), length(StartPoint:dim(RefDistMat)[1])), dimnames = ArrayDimNames)

    if (IgnorePrompts==FALSE){
      print(x = "Distance matrix rownames appear to be unique and have been used as IDs for DataDump file names or array dimensions.")
    }
  } else {
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(RefDistMat)))
    Prov.Array.Res <- array(NA, dim = c(length(Latways), length(Longways), length(StartPoint:dim(RefDistMat)[1])), dimnames = ArrayDimNames)

    if (IgnorePrompts==FALSE){
      print(x = "Distance matrix rownames appear to have duplicates, so order number has been used as IDs for DataDump file names or array dimensions.")
    }
  }



  IdentificationRange <- matrix(NA, nrow = dim(RefDistMat)[1], ncol = 3, dimnames = list(RefIDs, c('ID', '%RefDistribution', 'ProvencaningArea(m)')))
  for (dd in StartPoint:dim(RefDistMat)[1]){

    #dd <- 1

    if (PrintProg==TRUE){

      svMisc::progress(value = dd, max.value = dim(RefDistMat)[1], progress.bar = TRUE)
      #Sys.sleep(0.01)
      if (dd == dim(RefDistMat)[1]) cat("Done!\n")
    }

    #sum(is.infinite(RefDistMat))
    # extracting distances to this itteration's target specimen only
    Data.Dist <- RefDistMat[-dd,dd]


    #creating an empty object to be populated by results
    CoordsHeat <- NULL
    NarrowedRanges <- c(NA, NA)
    CorMatrixRes <- matrix(NA, nrow = length(Latways), ncol = length(Longways))
    rownames(CorMatrixRes) <- Latways
    colnames(CorMatrixRes) <- Longways
    #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
    for (I in 1:length(Longways)){
      for (J in 1:length(Latways)){
        i <- Longways[I]
        j <- Latways[J]


        coord <- c(j, i)

        GeographicDist <- GeoDist2Point(RefLatLongs = LatLongs[-dd,], TargetLatLong = coord)

        Data.Dist[is.infinite(Data.Dist)] <- 1

        if (Method=='Spearman'){
          #running the correlation to generate r
          CorRes <- suppressWarnings(stats::cor.test(x = Data.Dist, y = GeographicDist, method = "spearman"))
        } else if (Method=='Pearson'){
          CorRes <- suppressWarnings(stats::cor.test(x = Data.Dist, y = GeographicDist, method = "pearson"))
        }


        CorMatrixRes[J,I] <- CorRes$estimate


      }
    }


    if (sum(rownames(RefDistMat)==c(1:length(rownames(RefDistMat))))==length(rownames(RefDistMat))){
      SpecimenID <- dd
    } else if (length(unique(rownames(RefDistMat)))==length(rownames(RefDistMat))){
      SpecimenID <- rownames(RefDistMat)[dd]
    } else {
      SpecimenID <- dd
    }





    if (DataDump==TRUE){
      if (dd==StartPoint){
        BoundaryDDTot <- paste(DataDumpPath,"Individuals/", sep="")
        dir.create(BoundaryDDTot)
        utils::write.csv(CorMatrixRes, paste(BoundaryDDTot, SpecimenID, ".csv", sep=""))
      } else {
        utils::write.csv(CorMatrixRes, paste(BoundaryDDTot, SpecimenID, ".csv", sep=""))
      }

    }


    Prov.Array.Res[,,dd] <- CorMatrixRes


    #converting the results into a data frame
    CoordsHeat <- as.data.frame(cbind(expand.grid(rownames(CorMatrixRes), colnames(CorMatrixRes)), c(CorMatrixRes)), row.names = 1:length(c(CorMatrixRes)))



    #naming the variables
    names(CoordsHeat) <- c("Lats", "Longs", "Cor.Slope")

    CoordsHeatNum <- chr2nu(CoordsHeat$Cor.Slope)
    OriginLoc <- CoordsHeat[which(CoordsHeatNum==max(CoordsHeatNum)),]

    CoordsHeatscaled <- (CoordsHeatNum-min(CoordsHeatNum))/(max(CoordsHeatNum)-min(CoordsHeatNum))
    CoordsHeats <- grDevices::rgb(red = CoordsHeatscaled, green = 1-CoordsHeatscaled, blue = 1-CoordsHeatscaled)

    Select.95.conf <- which(chr2nu(CoordsHeat$Cor.Slope)>PlotValCor)


    if (length(Select.95.conf)!=0){

      ApproxOrigin <- CoordsHeat[Select.95.conf,]

      Latvar <- stats::var(chr2nu(ApproxOrigin$Lats))
      Longvar <- stats::var(chr2nu(ApproxOrigin$Longs))

      EstDist <- NULL
      if (is.na(Latvar)){
        if (is.null(RefIDs)){
          SpecimenMessage <- paste("For the", dd, "specimen", sep=" ")
        } else {
          SpecimenMessage <- paste("For the", RefIDs[dd], "specimen", sep=" ")
        }

        RthresMessage <- paste("A provenancing region was not identified at this r threshold. R threshold set to:", PlotValCor, sep = "   ")
        MaxCorMessage <- paste("The maximum correlation value acheived was:", max(chr2nu(CoordsHeat$Cor)), sep = "   ")
        WarningMessage <- "This may be because the resolution used is too low so the likely origin region has been overlooked or alternatively the specimen could not be successfully identified because the reference material does not suffieciently reflect the morphology of the unknown specimen. Please adjust the sampling resolution by changing the R.Samp argument or change the r threshold or consider the specimen unidentifiable to a given region."
        warning(paste(SpecimenMessage, RthresMessage, MaxCorMessage, WarningMessage, sep = " "))

      } else if (Latvar==0 || Longvar==0){
        if (is.null(RefIDs)){
          SpecimenMessage <- paste("For the", dd, "specimen", sep=" ")
        } else {
          SpecimenMessage <- paste("For the", RefIDs[dd], "specimen", sep=" ")
        }
        WarningMessage <- paste(SpecimenMessage, "The resolution used for identifying a region of identification is too low to plot a polygon of the likely region of origin. Therefore, the grid squares that were identified by the r threshold as a possible region of origin have been highlighted. Please set the R.Samp argument to a higher value if a polygon of the most likely region of origin is desired.")
        warning(WarningMessage)

        graphics::points(x = as.character(ApproxOrigin$Long)[1], y = as.character(ApproxOrigin$Lat)[1], pch=21, bg=transpar(Colour ='gold', alpha = dim(RefDistMat)[1]/10),  cex=dim(ApproxOrigin)[1])
        graphics::points(x = as.character(ApproxOrigin$Long), y = as.character(ApproxOrigin$Lat), pch=21, bg=transpar(Colour ='gold', alpha = dim(RefDistMat)[1]/10),  cex=dim(ApproxOrigin)[1])

      } else {
        contour.95 <-  Construct_contour(ApproxOrigin[,1:2])
        #polygon(contour.95$x, contour.95$y, col=transpar(Colour ='gold', alpha = dim(RefDistMat)[1]/10), border='black', lwd=1)

        EstimatedPol <-  sp::Polygon(cbind(contour.95$x, contour.95$y))
        EstimatedPols <-  sp::Polygons(list(EstimatedPol),1)
        EstSpatPol <-  sp::SpatialPolygons(list(EstimatedPols))
        #plot(DistSpatPol)
        #plot(EstSpatPol, add=TRUE)
        sp::proj4string(EstSpatPol) <- "+proj=longlat +ellps=WGS84"


        EstDist <- raster::intersect(EstSpatPol, DistSpatPol)

      }

      if (is.null(EstDist)){
        NarrowedRanges <- c(0, 0.00001) ###### check for better resolution, fails when resolution not high enough to put identification area within ref distribution area
      } else if (attr(class(EstDist), "package")!='sp'){
        NarrowedRanges <- c(0, 0.00001)
      } else {
        EstDist$area <- raster::area(EstDist)
        NarrowedRanges <- c(EstDist$area/DistSpatPol$area, EstDist$area)
      }

    } else {
      NarrowedRanges <- c(1, DistSpatPol$area)
    }

    IdentificationRange[dd,] <- c(as.character(RefIDs[dd]),NarrowedRanges)

    if (DataDump==TRUE){
      filename <- paste(DataDumpPath, "DataDump.csv", sep="")

      utils::write.csv(IdentificationRange, filename)
    }
  }


  Boundary.Array <- Prov.Array.Res

  for (A in 1:dim(Prov.Array.Res)[3]){
    Boundary.Array[,,A] <- BoundaryCalculation(MapMatrix = Prov.Array.Res[,,A], PlotValCor = PlotValCor)
  }


  Prov.matrix.Res <- matrix(NA, nrow = length(Latways), ncol = length(Longways))
  rownames(Prov.matrix.Res) <- Latways
  colnames(Prov.matrix.Res) <- Longways

  for (li in 1:length(Latways)){
    for (lj in 1:length(Longways)){
      Prov.matrix.Res[li,lj] <- sum(Boundary.Array[li,lj,])
    }
  }


  BoundaryCoordsResTable <- as.data.frame(cbind(expand.grid(rownames(Prov.matrix.Res), colnames(Prov.matrix.Res)), c(Prov.matrix.Res)), row.names = 1:length(c(CorMatrixRes)))

  BoundScaled <- (max(chr2nu(BoundaryCoordsResTable[,3]))-chr2nu(BoundaryCoordsResTable[,3]))/(max(chr2nu(BoundaryCoordsResTable[,3]))-min(chr2nu(BoundaryCoordsResTable[,3])))


  maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="red", bg=graphics::par(bg="white"))#, add=T, lwd=4)

  graphics::points(x = as.character(BoundaryCoordsResTable[,2]), y = as.character(BoundaryCoordsResTable[,1]), pch=15, col=grDevices::grey(BoundScaled),  cex=2)

  graphics::points(x = as.character(LatLongs$Long), y = as.character(LatLongs$Lat), pch=23, bg='blue',  cex=1)
  maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="red", bg=graphics::par(bg="white"), add=T, lwd=4)
  maps::map.axes()


  if (DataDump==TRUE){
    utils::write.csv(rownames(Prov.matrix.Res), paste(DataDumpPath, "LatKey.csv", sep = ""))
    utils::write.csv(colnames(Prov.matrix.Res), paste(DataDumpPath, "LongsKey.csv", sep = ""))

    return(IdentificationRange)
  } else {
    return(list(RawCorData=Prov.Array.Res, ProvenanceRange=IdentificationRange))

  }


}


#' Internal function for identifying boundaries from spatial correlation scores.
#'
#' This is an internal function to support the BoundaryFinder functions.
#' This function takes the correlation values calculated across the spatial grid and identifies which grid points are at the boundary of of the provenancing region.
#' These results then allow for the construction of compiled trait boundaries from spatial provenancing exercises.
#'
#' @return
#'
#' This function returns a matrix of logical TRUE/FALSE values identifying the boundary of the desired correlation value.
#'
#' @param MapMatrix a matrix of correlation values that correspond to grid reference locations
#' @inheritParams IDbyDistanceDistInput
#'
#' @keywords internal
#' @author Ardern Hulme-Beaman



BoundaryCalculation <- function(MapMatrix, PlotValCor){
  CorMatrixRes <- matrix(FALSE, nrow = dim(MapMatrix)[1], ncol = dim(MapMatrix)[2])
  rownames(CorMatrixRes) <- rownames(MapMatrix)
  colnames(CorMatrixRes) <- colnames(MapMatrix)

  Provencing.Mat <- MapMatrix>PlotValCor

  for (i in 1:dim(MapMatrix)[1]){
    for (j in 1:dim(MapMatrix)[2]){

      if (Provencing.Mat[i,j]==TRUE){
        if (i==1 & j==1){
          CorMatrixRes[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i+1,j]<3
        } else if (i==dim(MapMatrix)[1] & j==dim(MapMatrix)[2]){
          CorMatrixRes[i,j] <- Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i,j-1]<3
        } else if (i==dim(MapMatrix)[1] & j==1){
          CorMatrixRes[i,j] <- Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j+1]+Provencing.Mat[i,j+1]<3
        } else if (j==dim(MapMatrix)[2] & i==1){
          CorMatrixRes[i,j] <- Provencing.Mat[i+1,j]+Provencing.Mat[i+1,j-1]+Provencing.Mat[i,j-1]<3
        } else if (i==1){
          CorMatrixRes[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i+1,j]+Provencing.Mat[i,j+1]+Provencing.Mat[i,j-1]<5
        } else if (j==1){
          CorMatrixRes[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i+1,j]+Provencing.Mat[i-1,j+1]+Provencing.Mat[i-1,j]<5
        } else if (i==dim(MapMatrix)[1]){
          CorMatrixRes[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i-1,j+1]+Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i,j-1]<5
        } else if (j==dim(MapMatrix)[2]){
          CorMatrixRes[i,j] <- Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i,j-1]+Provencing.Mat[i+1,j-1]+Provencing.Mat[i+1,j]<5
        } else {
          CorMatrixRes[i,j] <- Provencing.Mat[i+1,j-1]+Provencing.Mat[i+1,j]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i,j-1]+Provencing.Mat[i,j+1]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j+1]<8
        }
      }


    }
  }


  return(CorMatrixRes)

}



#' Plots map with boundaries from data dump folder or R object
#'
#' This function plots the results taken from the \code{BoundaryFinder} function that have been dumped in
#' a specified folder or are stored in an array. These results are constructed from the correlation
#' values calculated across the spatial grid and identifies which gird points are at the boundary of
#' the provenancing region. These results then allows for the construction of compiled trait boundaries
#' from spatial provenancing exercises.
#' @param Path the directory where boundary data files have been dumped from the \code{BoundaryFinder} function. Default is set to NA.
#' @param DataDump if set to TRUE then the function uses files stored in the datadump folder as defined by \code{PATH}. If set to FALSE an array generated from the \code{BoundaryFinder} functions needs to be provided. Default is set to TRUE.
#' @param RawCorArray an array of correlation values for each reference specimen where columns and rows correspond with longitude and latitude and the third dimension represents each individual reference specimen. Default is set to NA.
#' @param MapLinesWd is a single numeric value to set the width of the map lines of landmass boundaries.
#' @param TileSize is a single numeric value to set the size of the pixels/tiles that will form the boundary.
#' @param RefPchSize is a single numeric value to set the size of the points marking where reference specimens are located.
#' @param plotLong is a vector of numeric values for all the Longitude values sampled in the boundary finding exercise. Unlike the LongRange argument of other functions where just the maximum and minimum values are specified, this vectors should contain all Longitude value expected. As a result, this numeric vector should match the column names of the data dump files.
#' @param plotLat is a vector of numeric values for all the Latitude values sampled in the boundary finding exercise. Unlike the LatRange argument of other functions where just the maximum and minimum values are specified, this vectors should contain all Latitude value expected. As a result, this numeric vector should match the row names of the data dump files.
#' @param MapExpansion is a vector of two numeric values that define an increase in plotting area of the map (but not the area for which the boundary calculations have been carried out). The first value is added to the  This therefore allows the
#' @param RefPchCol is the colour to be used for the reference specimen location points
#' @param BoundaryHue is a vector of 2 elements. The first element sets the colour hue value on a scale of 0 to 1 for a hsv function. The second value sets a transparancy level between 0 and 1, 0 being completely transparent.
#' @inheritParams BoundaryFinder
#' @inheritParams IDbyDistanceDistInput
#'
#' @details The map plotting of this function makes use of the functions of the \code{maps} package.
#'
#' @section Citations:
#'
#' Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
#' Enhancements by Thomas P Minka and Alex Deckmyn. (2017). maps: Draw Geographical Maps. R
#' package version 3.2.0. https://CRAN.R-project.org/package=maps
#'
#' @author Ardern Hulme-Beaman
#' @examples
#' FteydeaThres <- IDbyDistanceDistInputCCV(LatLongs = Fteydea$Info[,2:3],
#'                DistDataMat = Fteydea$SongDisMat,
#'                Verbose = TRUE,
#'                ProvConfidence = .95,
#'                PrintProg = FALSE)
#' PlotBoundaries(PlotValCor = FteydeaThres$`Provenancing.Correlation.95%.Confidence`,
#'                DataDump = FALSE,
#'                RawCorArray = Fteydea$Total.Boundary$RawCorData,
#'                plotLong = colnames(Fteydea$Total.Boundary$RawCorData),
#'                plotLat = rownames(Fteydea$Total.Boundary$RawCorData),
#'                LatLongs = Fteydea$Info[,2:3],
#'                TileSize = 2.5,
#'                BoundaryHue = c(0.5,1))
#'
#' @export


PlotBoundaries <- function(PlotValCor, DataDump=TRUE, Path=NA, RawCorArray=NA , MapLinesWd=1, TileSize=4, plotLong, plotLat, MapExpansion=c(0,0), LatLongs, RefPchCol='blue',  RefPchSize=1, BoundaryHue = c(1,1)){

  colnames(LatLongs) <- c('Lats', 'Longs')
  LatLongs <- as.data.frame(LatLongs)
  plotLong <- as.numeric(plotLong)
  plotLat <- as.numeric(plotLat)

  if (DataDump==TRUE){
    BoundaryCoords <- list.files(Path)

    #dimtest <- read.csv(paste(DataDumpPath, "Individuals/", BoundaryCoords[1], sep=""))[,-1]
    dimtest <- utils::read.csv(paste(Path, BoundaryCoords[1], sep=""))[,-1]

    Prov.Array.Res <- array(NA, dim = c(dim(dimtest)[1], dim(dimtest)[2], length(BoundaryCoords)))

    for (bb in 1:length(BoundaryCoords)){
      #bb <- 10
      ContentBoundary <- utils::read.csv(paste(Path,BoundaryCoords[bb], sep=""))[,-1]

      BoundaryMat <- as.matrix(ContentBoundary)
      Prov.Array.Res[,,bb] <- BoundaryCalculation(MapMatrix = BoundaryMat, PlotValCor = PlotValCor)
    }
  } else {
    Prov.Array.Res <- RawCorArray

    for (bb in 1:dim(Prov.Array.Res)[3]){

      BoundaryMat <- as.matrix(RawCorArray[,,bb])
      Prov.Array.Res[,,bb] <- BoundaryCalculation(MapMatrix = BoundaryMat, PlotValCor = PlotValCor)
    }
  }




  rownames(Prov.Array.Res) <- plotLat
  colnames(Prov.Array.Res) <- plotLong


  Prov.matrix.Res <- matrix(NA, nrow = length(plotLat), ncol = length(plotLong))
  rownames(Prov.matrix.Res) <- plotLat
  colnames(Prov.matrix.Res) <- plotLong

  for (li in 1:length(plotLat)){
    for (lj in 1:length(plotLong)){
      Prov.matrix.Res[li,lj] <- sum(Prov.Array.Res[li,lj,])
    }
  }



  Longways <- c(min(as.numeric(plotLong))-MapExpansion[2], max(as.numeric(plotLong))+MapExpansion[2])
  Latways <- c(min(as.numeric(plotLat))-MapExpansion[1], max(as.numeric(plotLat))+MapExpansion[1])



  BoundaryCoordsResTable <- as.data.frame(cbind(expand.grid(rownames(Prov.matrix.Res), colnames(Prov.matrix.Res)), c(Prov.matrix.Res)), row.names = 1:length(c(Prov.matrix.Res)))

  BoundScaled <- (max(chr2nu(BoundaryCoordsResTable[,3]))-chr2nu(BoundaryCoordsResTable[,3]))/(max(chr2nu(BoundaryCoordsResTable[,3]))-min(chr2nu(BoundaryCoordsResTable[,3])))


  maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="red", bg=graphics::par(bg="white"))#, add=T, lwd=4)

  BoundaryScale <- grDevices::hsv(h = BoundaryHue[1], v = 1, s = 1-BoundScaled, alpha = BoundaryHue[2])


  graphics::points(x = as.character(BoundaryCoordsResTable[,2]), y = as.character(BoundaryCoordsResTable[,1]), pch=15, col=BoundaryScale,  cex=TileSize)


  graphics::points(x = as.character(LatLongs$Longs), y = as.character(LatLongs$Lats), pch=23, bg=RefPchCol,  cex=RefPchSize)
  maps::map("world", xlim=c(min(plotLong), max(plotLong)), ylim=c(min(plotLat), max(plotLat)), interior=FALSE, col="grey45", bg=graphics::par(bg="white"), add=T, lwd=MapLinesWd)

  maps::map.axes()


}





#' Calculates trait boundary overlaps and statistics from Boundaryfinding R object
#'
#' This function takes the results of the \code{BoundaryFinder} function that have been
#' stored in an array and returns the statistics of the area coverage of boundaries
#' found for each specimen based on the provided r threshold. This function is particularly
#' useful when using the parallel processing function \code{BoundaryFinderPar}, which does not
#' return the Provenancing area values or percentages that are returned by the \code{BoundaryFInder}
#' function when the \code{verbose} aregument is set to \code{TRUE}.
#' @inheritParams PlotBoundaries
#' @inheritParams BoundaryFinder
#' @inheritParams IDbyDistanceDistInput
#'
#'
#' @return This function returns a list of two objects: 1. the area of the convex hull of the distribution of samples and 2. a dataframe of two columns with the area returned of reach specimen in the identification process and the corresponding percentage that area represents of the sample distribution.
#'
#'
#' @author Ardern Hulme-Beaman
#' @examples
#' FteydeaThres <- IDbyDistanceDistInputCCV(LatLongs = Fteydea$Info[,2:3],
#'                DistDataMat = Fteydea$SongDisMat,
#'                Verbose = TRUE,
#'                ProvConfidence = .95,
#'                PrintProg = FALSE)
#' TraitBoundaryStats(PlotValCor = FteydeaThres$`Provenancing.Correlation.95%.Confidence`,
#'                RawCorArray = Fteydea$Total.Boundary$RawCorData,
#'                LatLongs = Fteydea$Info[,2:3],
#'                RefIDs = Fteydea$Info[,1]
#'                )
#'
#' @export


TraitBoundaryStats <- function(RawCorArray, LatLongs, RefIDs, PlotValCor){

  SpecimenLoc <- cbind(chr2nu(LatLongs$Long), chr2nu(LatLongs$Lat))
  DistPol <-  sp::Polygon(SpecimenLoc[grDevices::chull(SpecimenLoc),])
  DistPols <-  sp::Polygons(list(DistPol),1)
  DistSpatPol <-  sp::SpatialPolygons(list(DistPols))
  sp::proj4string(DistSpatPol) <- "+proj=longlat +ellps=WGS84"
  DistSpatPol$area <- raster::area(DistSpatPol)

  ResTableDimNames <- list(RefIDs, c('Percent.Overlap', 'Area.Overlap.in.m'))
  ResTable <- matrix(data = NA,nrow = dim(RawCorArray)[3], ncol = 2, dimnames = ResTableDimNames)

  for (i in 1:dim(RawCorArray)[3]){
    #i <- 1

    CoordsHeat <- cbind(expand.grid(dimnames(RawCorArray)[[1]], dimnames(RawCorArray)[[2]]), c(RawCorArray[,,i]))

    colnames(CoordsHeat) <- c('Lats', 'Longs', 'Cor')

    Select.95.conf <- which(CoordsHeat$Cor>PlotValCor)


    if (length(Select.95.conf)!=0){

      ApproxOrigin <- CoordsHeat[Select.95.conf,]

      Latvar <- stats::var(chr2nu(ApproxOrigin$Lats))
      Longvar <- stats::var(chr2nu(ApproxOrigin$Longs))

      EstDist <- NULL
      if (is.na(Latvar)){
        if (is.null(RefIDs)){
          SpecimenMessage <- paste("For the", i, "specimen", sep=" ")
        } else {
          SpecimenMessage <- paste("For the", RefIDs[i], "specimen", sep=" ")
        }

        RthresMessage <- paste("A provenancing region was not identified at this r threshold. R threshold set to:", PlotValCor, sep = "   ")
        MaxCorMessage <- paste("The maximum correlation value acheived was:", max(chr2nu(CoordsHeat$Cor)), sep = "   ")
        WarningMessage <- "This may be because the resolution used is too low so the likely origin region has been overlooked or alternatively the specimen could not be successfully identified because the reference material does not suffieciently reflect the morphology of the unknown specimen. Please adjust the sampling resolution by changing the R.Samp argument or change the r threshold or consider the specimen unidentifiable to a given region."
        warning(paste(SpecimenMessage, RthresMessage, MaxCorMessage, WarningMessage, sep = " "))

      } else if (Latvar==0 || Longvar==0){
        if (is.null(RefIDs)){
          SpecimenMessage <- paste("For the", i, "specimen", sep=" ")
        } else {
          SpecimenMessage <- paste("For the", RefIDs[i], "specimen", sep=" ")
        }
        WarningMessage <- paste(SpecimenMessage, "The resolution used for identifying a region of identification is too low to plot a polygon of the likely region of origin. Therefore, the grid squares that were identified by the r threshold as a possible region of origin have been highlighted. Please set the R.Samp argument to a higher value if a polygon of the most likely region of origin is desired.")
        warning(WarningMessage)

      } else {
        contour.95 <-  Construct_contour(ApproxOrigin[,1:2])
        #polygon(contour.95$x, contour.95$y, col=transpar(Colour ='gold', alpha = 10), border='black', lwd=1)

        EstimatedPol <-  sp::Polygon(cbind(contour.95$x, contour.95$y))
        EstimatedPols <-  sp::Polygons(list(EstimatedPol),1)
        EstSpatPol <-  sp::SpatialPolygons(list(EstimatedPols))
        #plot(DistSpatPol)
        #plot(EstSpatPol, add=TRUE)
        sp::proj4string(EstSpatPol) <- "+proj=longlat +ellps=WGS84"


        EstDist <- raster::intersect(EstSpatPol, DistSpatPol)

      }

      if (is.null(EstDist)){
        NarrowedRanges <- c(0, 0) ###### check for better resolution, fails when resolution not high enough to put identification area within ref distribution area
      } else if (attr(class(EstDist), "package")!='sp'){
        NarrowedRanges <- c(0, 0)
      } else {
        EstDist$area <- raster::area(EstDist)
        NarrowedRanges <- c(EstDist$area/DistSpatPol$area, EstDist$area)
      }

    } else {
      NarrowedRanges <- c(1, DistSpatPol$area)
    }

    ResTable[i,] <- NarrowedRanges
  }

  ResTableDF <- as.data.frame(ResTable)

  ReturnRes <- list(Sample.Distribution.Convex.Hull.Area.m=DistSpatPol$area, Returned.Provenance.Range=ResTableDF)

  return(ReturnRes)
}





