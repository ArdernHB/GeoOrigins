


#' Returns geographic boundaries identified from trait distances
#'
#' This function takes the distances among all reference specimens and uses a process similar to
#' leave-one-out correct cross-validation to identify likely spatial trait boundaries. This method
#' only allows for the inclusion of specimens with specific known locations (i.e. specific
#' latitude-longitude coordinates). This functions requires a distance input, which allows the user
#' to input any desired dissimilarity or distance metrics as appropriate for the data.
#' @param Ref.Dist.mat is a square matrix of pairwise distances among all reference specimens
#' @param Expand.map is a vector of 2 elements for expanding the plotting region of the map. The first element expands the latitudinal area and the second element expands the longitudinal area.
#' @param DataDump determines whether a datadump is used, which is recommended for large datasets that may take time to computeand when the user might want to interupt the process mid way. The default is set to TRUE. If it's set to FALSE then the results are output as an array where rows and columns reflect latitude longitude and the third dimension are individual specimens in the reference dataset.
#' @param DataDumpPath is a directory path indicating where the data for constructing the boundaries can be dumped.
#' @param startpoint is an iteger that denotes the specimen to start the process on. As the method needs to cycle through all known specimens in the database, this can take some time. If the process is stopped for whatever reason the process can be picked up again by adjusting the startpoint to the specimen number that it had previously finished on.
#' @param Ref.IDs is a vector of the unique identifiers for each of the reference specimens in the reference dataset. These values will be used for naming the files in the datadump fololder and also for matching up with the returned summary results. Default is set to NULL and if this is not populated then reference data is worked through consecutively, naming the datadump files in consecutive order.
#' @param ignore.prompts default is set to FALSE, but if set to TRUE queries such as those confirming the location of the datadump will be surpressed.
#' @inheritParams IDbyDistance.DistInput
#' @return If verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If verbose is FALSE then only those grid references with the highest correlation values are returned.
#' @details The map plotting of this function makes use of the functions of the \code{maps} package.
#' Spatial polygons were constructed using functions from the \code{sp} package and calculations of polygon intersects and areas were carried out using functions from the \code{raster} package.
#' When the \code{print.prog} is set to TRUE, the \code{progress} function of the \code{svMisc} package is used.
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
#' @keywords Spatial provenancing
#' @keywords Spatial identification
#' @author Ardern Hulme-Beaman
#' @export



BoundaryFinder <- function(Lat.Longs, Ref.Dist.mat=matrix(), LongRange, LatRange, Range.Samp=10, print.prog=TRUE, plot.Val.cor, Expand.map=c(0,0), DataDump=TRUE, DataDumpPath=NA, startpoint=1, Ref.IDs=NULL, ignore.prompts=FALSE){
  #Lat.Longs = cbind(RatLatLongData$Lat, RatLatLongData$Long); Ref.Dist.mat = RatDistMat;LongRange = Long.Range; LatRange = Lat.Range; Range.Samp = R.Samp; plot.Val.cor = rThres$`Provenancing.Correlation.95%.Confidence`; Expand.map = c(0,0); Ref.IDs = RatLatLongData$ID; DataDump = FALSE
  #Lat.Longs = NEFtInfo[,2:3]; Ref.Dist.mat = NEDisMat; DataDump = FALSE; LongRange = NEFt.Long.Range; LatRange = NEFt.Lat.Range; Range.Samp = NEFtR.Samp; plot.Val.cor = NEFtrThres$`Provenancing.Correlation.95%.Confidence`;Ref.IDs = NEFtInfo$ID

  if (ignore.prompts==FALSE){
    if (DataDump==TRUE){
      if (is.na(DataDumpPath)){
        readline(prompt = "\n Datadump directory not set. Data will be dumped in current working directory. \n Press [Esc] now to exit function and change DataDumpPath argument or press enter to continue.")
      } else {
        if ((length(list.files(DataDumpPath))+length(list.dirs(DataDumpPath)))>1){
          QueryDD <- paste("\n Datadump directory",  DataDumpPath, " contains existing files. \n Press [Esc] now to exit function and change DataDumpPath argument or press enter to continue.")
          readline(prompt = QueryDD)
        } else {
          QueryDD <- paste("\n DataDumpPath directory set to",  DataDumpPath, "\n Press [Esc] now to exit function and change DataDumpPath argument or press enter to continue.")
          readline(prompt = QueryDD)
        }

      }
    } else {
      readline(prompt = "\n WARNING DataDump is set to FALSE. Results will be sorted in an array and if the process is stopped data will be lost. \n Press [Esc] now to exit function and change DataDumpPath argument or press enter to continue.")
    }
  }



  #making Lat.Longs a dataframe
  Lat.Longs <- as.data.frame(Lat.Longs)
  colnames(Lat.Longs) <- c("Lats", "Longs")

  #creating a range that will cover the whole geographic area of interest
  #this is for the function to loop through later on
  if (length(Range.Samp)==1){
    Long.Samp <- Range.Samp
    Lat.Samp <- Range.Samp
  } else if (length(Range.Samp)==2){
    Long.Samp <- Range.Samp[2]
    Lat.Samp <- Range.Samp[1]
  } else if (length(Range.Samp)>2){
    stop("too many dimensions in Range.Samp")
  }

  Long.Range.Steps <- (LongRange[2]-LongRange[1])/(Long.Samp-2)
  Lat.Range.Steps <- (LatRange[2]-LatRange[1])/(Lat.Samp-2)

  #this output for Lat/Long ways provides what the loop should sequence through
  Longways <- c(LongRange[1]-Long.Range.Steps, seq(LongRange[1], LongRange[2], by = Long.Range.Steps), LongRange[2]+Long.Range.Steps)
  Latways <- c(LatRange[1]-Lat.Range.Steps, seq(LatRange[1], LatRange[2], by = Lat.Range.Steps), LatRange[2]+Lat.Range.Steps)

  #plotting a map that will later be populated with polygons for each jackknifing iteration in the dd loop below
  if (sum(Expand.map>0)>0){
    maps::map("world", xlim=c(min(Longways)-Expand.map[2], max(Longways)+Expand.map[2]), ylim=c(min(Latways)-Expand.map[1], max(Latways)+Expand.map[1]), interior=FALSE, col="black", bg=graphics::par(bg="white"))
    graphics::points(x = as.character(Lat.Longs$Long), y = as.character(Lat.Longs$Lat), pch=23, bg='blue',  cex=1)
  } else {
    maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"))
    graphics::points(x = as.character(Lat.Longs$Long), y = as.character(Lat.Longs$Lat), pch=23, bg='blue',  cex=1)
  }

  SpecimenLoc <- cbind(chr2nu(Lat.Longs$Long), chr2nu(Lat.Longs$Lat))
  DistPol <-  sp::Polygon(SpecimenLoc[grDevices::chull(SpecimenLoc),])
  DistPols <-  sp::Polygons(list(DistPol),1)
  DistSpatPol <-  sp::SpatialPolygons(list(DistPols))
  sp::proj4string(DistSpatPol) <- "+proj=longlat +ellps=WGS84"
  DistSpatPol$area <- raster::area(DistSpatPol)


  if (is.null(Ref.IDs)){
    Ref.IDs <- 1:dim(Ref.Dist.mat)[1]
  }

  if (sum(rownames(Ref.Dist.mat)==c(1:length(rownames(Ref.Dist.mat))))==length(rownames(Ref.Dist.mat))){
    ArrayDimNames <- list(Latways, Longways, rownames(Ref.Dist.mat))
    Prov.Array.Res <- array(NA, dim = c(length(Latways), length(Longways), length(startpoint:dim(Ref.Dist.mat)[1])), dimnames = ArrayDimNames)

    if (ignore.prompts==FALSE){
      readline(prompt = "\n Distance matrix rownames appear match numeric order, \n so order number has been used as IDs for DataDump file names or array dimensions. \n \n Press enter to continue.")
    }
  } else if (length(unique(rownames(Ref.Dist.mat)))==length(rownames(Ref.Dist.mat))){
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(Ref.Dist.mat)))
    Prov.Array.Res <- array(NA, dim = c(length(Latways), length(Longways), length(startpoint:dim(Ref.Dist.mat)[1])), dimnames = ArrayDimNames)

    if (ignore.prompts==FALSE){
      readline(prompt = "\n Distance matrix rownames appear to be unique \n and have been used as IDs for DataDump file names or array dimensions. \n \n Press enter to continue.")
    }
  } else {
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(Ref.Dist.mat)))
    Prov.Array.Res <- array(NA, dim = c(length(Latways), length(Longways), length(startpoint:dim(Ref.Dist.mat)[1])), dimnames = ArrayDimNames)

    if (ignore.prompts==FALSE){
      readline(prompt = "\n Distance matrix rownames appear to have duplicates, \n so order number has been used as IDs for DataDump file names or array dimensions. \n \n Press enter to continue.")
    }
  }



  IdentificationRange <- matrix(NA, nrow = dim(Ref.Dist.mat)[1], ncol = 3, dimnames = list(Ref.IDs, c('ID', '%RefDistribution', 'ProvencaningArea(m)')))
  for (dd in startpoint:dim(Ref.Dist.mat)[1]){

    #dd <- 1

    if (print.prog==TRUE){

      svMisc::progress(value = dd, max.value = dim(Ref.Dist.mat)[1], progress.bar = TRUE)
      Sys.sleep(0.01)
      if (dd == dim(Ref.Dist.mat)[1]) cat("Done!\n")
    }

    #sum(is.infinite(Ref.Dist.mat))
    # extracting distances to this itteration's target specimen only
    Data.Dist <- Ref.Dist.mat[-dd,dd]


    #creating an empty object to be populated by results
    CoordsHeat <- NULL
    NarrowedRanges <- c(NA, NA)
    cor.matrix.Res <- matrix(NA, nrow = length(Latways), ncol = length(Longways))
    rownames(cor.matrix.Res) <- Latways
    colnames(cor.matrix.Res) <- Longways
    #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
    for (I in 1:length(Longways)){
      for (J in 1:length(Latways)){
        i <- Longways[I]
        j <- Latways[J]


        coord <- c(j, i)

        Geographic.Dists <- Geo.Dist2Point(RefLatLongs = Lat.Longs[-dd,], TargetLatLong = coord)

        Data.Dist[is.infinite(Data.Dist)] <- 1

        #running the correlation to generate r
        CorRes <- stats::cor.test(x = Data.Dist, y = Geographic.Dists)

        cor.matrix.Res[J,I] <- CorRes$estimate


      }
    }


    if (sum(rownames(Ref.Dist.mat)==c(1:length(rownames(Ref.Dist.mat))))==length(rownames(Ref.Dist.mat))){
      SpecimenID <- dd
    } else if (length(unique(rownames(Ref.Dist.mat)))==length(rownames(Ref.Dist.mat))){
      SpecimenID <- rownames(Ref.Dist.mat)[dd]
    } else {
      SpecimenID <- dd
    }





    if (DataDump==TRUE){
      if (dd==startpoint){
        BoundaryDDTot <- paste(DataDumpPath,"Individuals/", sep="")
        dir.create(BoundaryDDTot)
        utils::write.csv(cor.matrix.Res, paste(BoundaryDDTot, SpecimenID, ".csv", sep=""))
      } else {
        utils::write.csv(cor.matrix.Res, paste(BoundaryDDTot, SpecimenID, ".csv", sep=""))
      }

    }


    Prov.Array.Res[,,dd] <- cor.matrix.Res


    #converting the results into a data frame
    CoordsHeat <- as.data.frame(cbind(expand.grid(rownames(cor.matrix.Res), colnames(cor.matrix.Res)), c(cor.matrix.Res)), row.names = 1:length(c(cor.matrix.Res)))



    #naming the variables
    names(CoordsHeat) <- c("Lats", "Longs", "Cor.Slope")

    CoordsHeatNum <- chr2nu(CoordsHeat$Cor.Slope)
    OriginLoc <- CoordsHeat[which(CoordsHeatNum==max(CoordsHeatNum)),]

    CoordsHeatscaled <- (CoordsHeatNum-min(CoordsHeatNum))/(max(CoordsHeatNum)-min(CoordsHeatNum))
    CoordsHeats <- grDevices::rgb(red = CoordsHeatscaled, green = 1-CoordsHeatscaled, blue = 1-CoordsHeatscaled)

    Select.95.conf <- which(chr2nu(CoordsHeat$Cor.Slope)>plot.Val.cor)


    if (length(Select.95.conf)!=0){

      ApproxOrigin <- CoordsHeat[Select.95.conf,]

      Latvar <- stats::var(chr2nu(ApproxOrigin$Lats))
      Longvar <- stats::var(chr2nu(ApproxOrigin$Longs))

      EstDist <- NULL
      if (is.na(Latvar)){
        if (is.null(Ref.IDs)){
          SpecimenMessage <- paste("For the", dd, "specimen", sep=" ")
        } else {
          SpecimenMessage <- paste("For the", Ref.IDs[dd], "specimen", sep=" ")
        }

        RthresMessage <- paste("A provenancing region was not identified at this r threshold. R threshold set to:", plot.Val.cor, sep = "   ")
        MaxCorMessage <- paste("The maximum correlation value acheived was:", max(chr2nu(CoordsHeat$Cor)), sep = "   ")
        WarningMessage <- "This may be because the resolution used is too low so the likely origin region has been overlooked or alternatively the specimen could not be successfully identified because the reference material does not suffieciently reflect the morphology of the unknown specimen. Please adjust the sampling resolution by changing the R.Samp argument or change the r threshold or consider the specimen unidentifiable to a given region."
        warning(paste(SpecimenMessage, RthresMessage, MaxCorMessage, WarningMessage, sep = " "))

      } else if (Latvar==0 || Longvar==0){
        if (is.null(Ref.IDs)){
          SpecimenMessage <- paste("For the", dd, "specimen", sep=" ")
        } else {
          SpecimenMessage <- paste("For the", Ref.IDs[dd], "specimen", sep=" ")
        }
        WarningMessage <- paste(SpecimenMessage, "The resolution used for identifying a region of identification is too low to plot a polygon of the likely region of origin. Therefore, the grid squares that were identified by the r threshold as a possible region of origin have been highlighted. Please set the R.Samp argument to a higher value if a polygon of the most likely region of origin is desired.")
        warning(WarningMessage)

        graphics::points(x = as.character(ApproxOrigin$Long)[1], y = as.character(ApproxOrigin$Lat)[1], pch=21, bg=transpar(Colour ='gold', alpha = dim(Ref.Dist.mat)[1]/10),  cex=dim(ApproxOrigin)[1])
        graphics::points(x = as.character(ApproxOrigin$Long), y = as.character(ApproxOrigin$Lat), pch=21, bg=transpar(Colour ='gold', alpha = dim(Ref.Dist.mat)[1]/10),  cex=dim(ApproxOrigin)[1])

      } else {
        contour.95 <-  Construct_contour(ApproxOrigin[,1:2])
        #polygon(contour.95$x, contour.95$y, col=transpar(Colour ='gold', alpha = dim(Ref.Dist.mat)[1]/10), border='black', lwd=1)

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

    IdentificationRange[dd,] <- c(as.character(Ref.IDs[dd]),NarrowedRanges)

    if (DataDump==TRUE){
      filename <- paste(DataDumpPath, "DataDump.csv", sep="")

      utils::write.csv(IdentificationRange, filename)
    }
  }


  Boundary.Array <- Prov.Array.Res

  for (A in 1:dim(Prov.Array.Res)[3]){
    Boundary.Array[,,A] <- BoundaryCalculation(Map.matrix = Prov.Array.Res[,,A], plot.Val.cor = plot.Val.cor)
  }


  Prov.matrix.Res <- matrix(NA, nrow = length(Latways), ncol = length(Longways))
  rownames(Prov.matrix.Res) <- Latways
  colnames(Prov.matrix.Res) <- Longways

  for (li in 1:length(Latways)){
    for (lj in 1:length(Longways)){
      Prov.matrix.Res[li,lj] <- sum(Boundary.Array[li,lj,])
    }
  }


  BoundaryCoordsResTable <- as.data.frame(cbind(expand.grid(rownames(Prov.matrix.Res), colnames(Prov.matrix.Res)), c(Prov.matrix.Res)), row.names = 1:length(c(cor.matrix.Res)))

  BoundScaled <- (max(chr2nu(BoundaryCoordsResTable[,3]))-chr2nu(BoundaryCoordsResTable[,3]))/(max(chr2nu(BoundaryCoordsResTable[,3]))-min(chr2nu(BoundaryCoordsResTable[,3])))


  maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="red", bg=graphics::par(bg="white"))#, add=T, lwd=4)

  graphics::points(x = as.character(BoundaryCoordsResTable[,2]), y = as.character(BoundaryCoordsResTable[,1]), pch=15, col=grDevices::grey(BoundScaled),  cex=2)

  graphics::points(x = as.character(Lat.Longs$Long), y = as.character(Lat.Longs$Lat), pch=23, bg='blue',  cex=1)
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


#' Internal function for identifies boundaries from spatial correlation scores.
#'
#' This function takes the correlation values calculated across the spatial grid and identifies which gird points are at the boundary of of the provenancing region.
#' These results then allows for the construction of compiled trait boundaries from spatial provenancing excercises.
#'
#' @return
#'
#' This function returns a matrix of logical TRUE/FALSE values identifying the boundary of the desired correlation value.
#'
#' @param Map.matrix is a matrix of correlation values that correspond to grid reference locations
#' @inheritParams IDbyDistance.DistInput
#'
#' @keywords Spatial provenancing; Boundary finder; internal.
#' @author Ardern Hulme-Beaman


BoundaryCalculation <- function(Map.matrix, plot.Val.cor){
  cor.matrix.Res <- matrix(FALSE, nrow = dim(Map.matrix)[1], ncol = dim(Map.matrix)[2])
  rownames(cor.matrix.Res) <- rownames(Map.matrix)
  colnames(cor.matrix.Res) <- colnames(Map.matrix)

  Provencing.Mat <- Map.matrix>plot.Val.cor

  for (i in 1:dim(Map.matrix)[1]){
    for (j in 1:dim(Map.matrix)[2]){

      if (Provencing.Mat[i,j]==TRUE){
        if (i==1 & j==1){
          cor.matrix.Res[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i+1,j]<3
        } else if (i==dim(Map.matrix)[1] & j==dim(Map.matrix)[2]){
          cor.matrix.Res[i,j] <- Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i,j-1]<3
        } else if (i==dim(Map.matrix)[1] & j==1){
          cor.matrix.Res[i,j] <- Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j+1]+Provencing.Mat[i,j+1]<3
        } else if (j==dim(Map.matrix)[2] & i==1){
          cor.matrix.Res[i,j] <- Provencing.Mat[i+1,j]+Provencing.Mat[i+1,j-1]+Provencing.Mat[i,j-1]<3
        } else if (i==1){
          cor.matrix.Res[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i+1,j]+Provencing.Mat[i,j+1]+Provencing.Mat[i,j-1]<5
        } else if (j==1){
          cor.matrix.Res[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i+1,j]+Provencing.Mat[i-1,j+1]+Provencing.Mat[i-1,j]<5
        } else if (i==dim(Map.matrix)[1]){
          cor.matrix.Res[i,j] <- Provencing.Mat[i,j+1]+Provencing.Mat[i-1,j+1]+Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i,j-1]<5
        } else if (j==dim(Map.matrix)[2]){
          cor.matrix.Res[i,j] <- Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i,j-1]+Provencing.Mat[i+1,j-1]+Provencing.Mat[i+1,j]<5
        } else {
          cor.matrix.Res[i,j] <- Provencing.Mat[i+1,j-1]+Provencing.Mat[i+1,j]+Provencing.Mat[i+1,j+1]+Provencing.Mat[i,j-1]+Provencing.Mat[i,j+1]+Provencing.Mat[i-1,j-1]+Provencing.Mat[i-1,j]+Provencing.Mat[i-1,j+1]<8
        }
      }


    }
  }


  return(cor.matrix.Res)

}



#' Plots map with boundaries from data dump folder
#'
#' This function plots the results of takes the BoundaryFinder function that have been dumped in
#' a specified folder or are stored in an array. These results are constructed from the correlation
#' values calculated across the spatial grid and identifies which gird points are at the boundary of
#' the provenancing region. These results then allows for the construction of compiled trait boundaries
#' from spatial provenancing excercises.
#' @param Path the directory where boundary data files have been dumped from the BoundaryFinder function. Default is set to NA.
#' @param DataDump if set to TRUE then the function uses files stored in the datadump folder as defined by \code{PATH}. If set to FALSE an array generated from the \code{BoundaryFinder} functions needs to be provided. Default is set to TRUE.
#' @param RawCorArray an array of correlation values for each reference specimen where columns and rows correspond with longitude and latitude and the third dimension represents each individual reference specimen. Default is set to NA.
#' @param MapLinesWd is a single numeric value to set the width of the map lines of landmass boundaries.
#' @param TileSize is a single numeric value to set the size of the pixels/tiles that will form the boundary.
#' @param Ref.Pch.Size is a single numeric value to set the size of the points marking where reference specimens are located.
#' @param plotLong is a vector of numeric values for all the Longitude values sampled in the boundary finding exercise. Unlike the LongRange argument of other functions where just the maximum and minimum values are specified, this vectors should contain all Longitude value expected. As a result, this numeric vector should match the column names of the data dump files.
#' @param plotLat is a vector of numeric values for all the Latitude values sampled in the boundary finding exercise. Unlike the LatRange argument of other functions where just the maximum and minimum values are specified, this vectors should contain all Latitude value expected. As a result, this numeric vector should match the row names of the data dump files.
#' @param MapExpansion is a vector of two numeric values that define an increase in plotting area of the map (but not the area for which the boundary calculations have been carried out). The first value is added to the  This therefore allows the
#' @param Ref.PtCol is the colour to be used for the reference specimen location points
#' @param BoundaryHue is a vector of 2 elements. The first element sets the colour hue value on a scale of 0 to 1 for a hsv function. The second value sets a transparancy level between 0 and 1, 0 being completely transparent.
#' @inheritParams BoundaryFinder
#' @inheritParams IDbyDistance.DistInput
#' @return matrix of logical TRUE/FALSE values identifying the boundary of the desired correlation value
#'
#' @details The map plotting of this function makes use of the functions of the \code{maps} package.
#'
#' @section Citations:
#'
#' Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
#' Enhancements by Thomas P Minka and Alex Deckmyn. (2017). maps: Draw Geographical Maps. R
#' package version 3.2.0. https://CRAN.R-project.org/package=maps
#'
#' @keywords Spatial provenancing
#' @keywords Boundary finder
#' @keywords Boundary plotting
#' @export


PlotBoundaries <- function(plot.Val.cor, DataDump=TRUE, Path=NA, RawCorArray=NA , MapLinesWd=1, TileSize=4, plotLong, plotLat, MapExpansion=c(0,0), Lat.Longs, Ref.PtCol='blue',  Ref.Pch.Size=1, BoundaryHue = c(1,1)){
  #plot.Val.cor = rThres$`Provenancing.Correlation.95%.Confidence`; RawCorArray = Boundaryfinding$RawCorData; plotLong = colnames(Boundaryfinding$RawCorData); plotLat = rownames(Boundaryfinding$RawCorData); Lat.Longs= cbind(RatLatLongData$Lat, RatLatLongData$Long)
  #DataDump=FALSE

  colnames(Lat.Longs) <- c('Lats', 'Longs')
  Lat.Longs <- as.data.frame(Lat.Longs)
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
      Prov.Array.Res[,,bb] <- BoundaryCalculation(Map.matrix = BoundaryMat, plot.Val.cor = plot.Val.cor)
    }
  } else {
    Prov.Array.Res <- RawCorArray

    for (bb in 1:dim(Prov.Array.Res)[3]){

      BoundaryMat <- as.matrix(RawCorArray[,,bb])
      Prov.Array.Res[,,bb] <- BoundaryCalculation(Map.matrix = BoundaryMat, plot.Val.cor = plot.Val.cor)
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


  graphics::points(x = as.character(Lat.Longs$Longs), y = as.character(Lat.Longs$Lats), pch=23, bg=Ref.PtCol,  cex=Ref.Pch.Size)
  maps::map("world", xlim=c(min(plotLong), max(plotLong)), ylim=c(min(plotLat), max(plotLat)), interior=FALSE, col="grey45", bg=graphics::par(bg="white"), add=T, lwd=MapLinesWd)

  maps::map.axes()


}






