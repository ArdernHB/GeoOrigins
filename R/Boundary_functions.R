


#' Returns geographic boundaries identified from trait distances
#'
#' This function takes the distances among all reference specimens and uses a process similar to
#' leave-one-out correct cross-validation to identify likely spatial trait boundaries. This method
#' only allows for the inclusion of specimens with specific known locations (i.e. specific
#' latitude-longitude coordinates). This functions requires a distance input, which allows the user
#' to input any desired dissimilarity or distance metrics as appropriate for the data.
#' @param Ref.Dist.mat is a square matrix of pairwise distances among all reference specimens
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



BoundaryFinder <- function(Lat.Longs=data.frame(Lats=c(), Longs=c()), Ref.Dist.mat=matrix(), LongRange=c(0,0), LatRange=c(0,0), Range.Samp=10, print.prog=TRUE, plot.Val.cor=as.numeric(X), Expand.map=c(Lat, Long), datadump=NA, startpoint=1){

  #Lat.Longs=Heat.Info; Ref.Dist.mat=HeatDistMat; LongRange = Long.Range; LatRange = Lat.Range; Range.Samp = R.Samp; print.prog=TRUE; plot.Val.cor=.66; Expand.map = c(.05, .02); datadump = "DataDump_Temp_folder_test/"; startpoint=1

  if (is.na(datadump)){
    readline(prompt = "\n Datadump directory not set. Data will be dumped in current working directory. \n Press [Esc] now to exit function and change datadump argument or press any key to continue.")
  } else {
    if ((length(list.files(datadump))+length(list.dirs(datadump)))>1){
      QueryDD <- paste("\n Datadump directory",  datadump, " contains existing files. \n Press [Esc] now to exit function and change datadump argument or press any key to continue.")
      readline(prompt = QueryDD)
    } else {
      QueryDD <- paste("\n Datadump directory set to",  datadump, "\n Press [Esc] now to exit function and change datadump argument or press any key to continue.")
      readline(prompt = QueryDD)
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
    maps::map("world", xlim=c(min(Longways)-Expand.map[2], max(Longways)+Expand.map[2]), ylim=c(min(Latways)-Expand.map[1], max(Latways)+Expand.map[1]), interior=FALSE, col="black", bg=par(bg="white"))
    points(x = as.character(Lat.Longs$Long), y = as.character(Lat.Longs$Lat), pch=23, bg='blue',  cex=1)
  } else {
    maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=par(bg="white"))
    points(x = as.character(Lat.Longs$Long), y = as.character(Lat.Longs$Lat), pch=23, bg='blue',  cex=1)
  }

  SpecimenLoc <- cbind(chr2nu(Lat.Longs$Long), chr2nu(Lat.Longs$Lat))
  DistPol <-  sp::Polygon(SpecimenLoc[chull(SpecimenLoc),])
  DistPols <-  sp::Polygons(list(DistPol),1)
  DistSpatPol <-  sp::SpatialPolygons(list(DistPols))
  DistSpatPol$area <- raster::area(DistSpatPol)

  Prov.Array.Res <- array(NA, dim = c(length(Latways), length(Longways), length(startpoint:dim(Ref.Dist.mat)[1])))

  IdentificationRange <- NULL
  for (dd in startpoint:dim(Ref.Dist.mat)[1]){

    #dd <- 3

    if (print.prog==TRUE){


      svMisc::progress(dd, progress.bar = TRUE)
      Sys.sleep(0.01)
      if (i == dim(Ref.Dist.mat)[1]) cat("Done!\n")
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


    if (sum(rownames(HeatDistMat)==c(1:length(rownames(HeatDistMat))))==length(rownames(HeatDistMat))){
      SpecimenID <- dd
      if (dd==1){
        print(" Distance matrix rownames appear match numeric order, \n so order number has been used as IDs for datadump file names")
      }
    } else if (length(unique(rownames(HeatDistMat)))==length(rownames(HeatDistMat))){
      SpecimenID <- rownames(HeatDistMat)[dd]
      if (dd==1){
        print(" Distance matrix rownames appear to be unique \n and have been used as IDs for datadump file names")
      }
    } else {
      SpecimenID <- dd
      if (dd==1){
        print(" Distance matrix rownames appear to have duplicates, \n so order number has been used as IDs for datadump file names")
      }
    }


    if (dd==startpoint){
      BoundaryDDTot <- paste(datadump,"Individuals/", sep="")
      dir.create(BoundaryDDTot)
      write.csv(cor.matrix.Res, paste(BoundaryDDTot, SpecimenID, ".csv", sep=""))
    } else {
      write.csv(cor.matrix.Res, paste(BoundaryDDTot, SpecimenID, ".csv", sep=""))
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

      EstDist <- 0
      if (dim(ApproxOrigin)[1]<5){

        points(x = as.character(ApproxOrigin$Long)[1], y = as.character(ApproxOrigin$Lat)[1], pch=21, bg=transpar(Colour ='gold', alpha = dim(Ref.Dist.mat)[1]/10),  cex=dim(ApproxOrigin)[1])
        points(x = as.character(ApproxOrigin$Long), y = as.character(ApproxOrigin$Lat), pch=21, bg=transpar(Colour ='gold', alpha = dim(Ref.Dist.mat)[1]/10),  cex=dim(ApproxOrigin)[1])

      } else {
        contour.95 <-  Construct_contour(ApproxOrigin[,1:2])
        #polygon(contour.95$x, contour.95$y, col=transpar(Colour ='gold', alpha = dim(Ref.Dist.mat)[1]/10), border='black', lwd=1)

        EstimatedPol <-  sp::Polygon(cbind(contour.95$x, contour.95$y))
        EstimatedPols <-  sp::Polygons(list(EstimatedPol),1)
        EstSpatPol <-  sp::SpatialPolygons(list(EstimatedPols))
        #plot(DistSpatPol)
        #plot(EstSpatPol, add=TRUE)


        EstDist <- raster::intersect(EstSpatPol, DistSpatPol)

      }

      if (is.null(EstDist)){
        NarrowedRanges <- c(0, 0.00001) ###### check for better solution, fails when resolution not high enough to put identification area within ref distribution area
      } else if (attr(class(EstDist), "package")!='sp'){
        NarrowedRanges <- c(0, 0.00001)
      } else {
        EstDist$area <- raster::area(EstDist)
        NarrowedRanges <- c(EstDist$area/DistSpatPol$area, EstDist$area)
      }

    } else {
      NarrowedRanges <- c(1, DistSpatPol$area)
    }

    IdentificationRange <- rbind(IdentificationRange, c(dd,NarrowedRanges))

    if (is.na(datadump)==FALSE){
      filename <- paste(datadump, "DataDump.csv", sep="")

      write.csv(IdentificationRange, filename)
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


  maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="red", bg=par(bg="white"))#, add=T, lwd=4)

  points(x = as.character(BoundaryCoordsResTable[,2]), y = as.character(BoundaryCoordsResTable[,1]), pch=15, col=grey(BoundScaled),  cex=2)

  points(x = as.character(Lat.Longs$Long), y = as.character(Lat.Longs$Lat), pch=23, bg='blue',  cex=1)
  maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="red", bg=par(bg="white"), add=T, lwd=4)
  maps::map.axes()

  write.csv(rownames(Prov.matrix.Res), paste(datadump, "LatKey.csv", sep = ""))
  write.csv(colnames(Prov.matrix.Res), paste(datadump, "LongsKey.csv", sep = ""))
  return(IdentificationRange)
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
#' a specified folder. These results are constructed from the correlation values calculated across
#' the spatial grid and identifies which gird points are at the boundary of of the provenancing
#' region. These results then allows for the construction of compiled trait boundaries from spatial
#' provenancing excercises.
#' @param Path the directory where boundary data files have been dumped from the BoundaryFinder function.
#' @param MapLinesWd is a single numeric value to set the width of the map lines of landmass boundaries.
#' @param TileSize is a single numeric value to set the size of the pixels/tiles that will form the boundary.
#' @param Ref.Pch.Size is a single numeric value to set the size of the points marking where reference specimens are located.
#' @param plotLong is a vector of numeric values for all the Longitude values sampled in the boundary finding exercise. Unlike the LongRange argument of other functions where just the maximum and minimum values are specified, this vectors should contain all Longitude value expected. As a result, this numeric vector should match the column names of the data dump files.
#' @param plotLat is a vector of numeric values for all the Latitude values sampled in the boundary finding exercise. Unlike the LatRange argument of other functions where just the maximum and minimum values are specified, this vectors should contain all Latitude value expected. As a result, this numeric vector should match the row names of the data dump files.
#' @param MapExpansion is a vector of two numeric values that define an increase in plotting area of the map (but not the area for which the boundary calculations have been carried out). The first value is added to the  This therefore allows the
#' @param Ref.PtCol is the colour to be used for the reference specimen location points
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


PlotBoundaries <- function(plot.Val.cor, Path, MapLinesWd=1, TileSize=4, plotLong, plotLat, MapExpansion=c(0,0), Lat.Longs=data.frame(Lats=c(), Longs=c()), Ref.PtCol='blue',  Ref.Pch.Size=1, BoundaryHue = c(1,1)){

  #CorVal=plot.Val.cor; plotLong= Long.Range; plotLat=Lat.Range; Path="NewBoundaryTrial_fresh/Individuals/"; MapLinesWd = 1; TileSize = 2; MapExpansion = c(.2,.1); Spec.Lat.Longs = Heat.Info

  BoundaryCoords <- list.files(Path)

  #dimtest <- read.csv(paste(datadump, "Individuals/", BoundaryCoords[1], sep=""))[,-1]
  dimtest <- read.csv(paste(Path, BoundaryCoords[1], sep=""))[,-1]

  Prov.Array.Res <- array(NA, dim = c(dim(dimtest)[1], dim(dimtest)[2], length(BoundaryCoords)))

  for (bb in 1:length(BoundaryCoords)){
    #bb <- 10
    ContentBoundary <- read.csv(paste(Path,BoundaryCoords[bb], sep=""))[,-1]

    BoundaryMat <- as.matrix(ContentBoundary)
    Prov.Array.Res[,,bb] <- BoundaryCalculation(Map.matrix = BoundaryMat, plot.Val.cor = plot.Val.cor)
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



  Longways <- c(min(plotLong)-MapExpansion[2], max(plotLong)+MapExpansion[2])
  Latways <- c(min(plotLat)-MapExpansion[1], max(plotLat)+MapExpansion[1])



  BoundaryCoordsResTable <- as.data.frame(cbind(expand.grid(rownames(Prov.matrix.Res), colnames(Prov.matrix.Res)), c(Prov.matrix.Res)), row.names = 1:length(c(Prov.matrix.Res)))

  BoundScaled <- (max(chr2nu(BoundaryCoordsResTable[,3]))-chr2nu(BoundaryCoordsResTable[,3]))/(max(chr2nu(BoundaryCoordsResTable[,3]))-min(chr2nu(BoundaryCoordsResTable[,3])))


  maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="red", bg=par(bg="white"))#, add=T, lwd=4)

  BoundaryScale <- hsv(h = BoundaryHue[1], v = 1, s = 1-BoundScaled, alpha = BoundaryHue[2])


  points(x = as.character(BoundaryCoordsResTable[,2]), y = as.character(BoundaryCoordsResTable[,1]), pch=15, col=BoundaryScale,  cex=TileSize)


  points(x = as.character(plot.Val.Lat.Longs$Longs), y = as.character(plot.Val.Lat.Longs$Lats), pch=23, bg=Ref.PtCol,  cex=Ref.Pch.Size)
  maps::map("world", xlim=c(min(plotLong), max(plotLong)), ylim=c(min(plotLat), max(plotLat)), interior=FALSE, col="grey45", bg=par(bg="white"), add=T, lwd=MapLinesWd)

  maps::map.axes()


}






