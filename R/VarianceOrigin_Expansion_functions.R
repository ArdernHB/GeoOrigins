



#' Heat mapping function to identify region of highest value
#'
#' This function takes a metric for spatial specimens or populations with defined locations
#' and calculates the region of highest values. This can be used with population variance
#' scores to identify regions of highest variance.
#' @param VarData a numeric vector of values (e.g. variance metrics) in corresponding order with the LatLongs.
#' @inheritParams IDbyDistanceRawData
#' @return If Verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If Verbose is FALSE then only those grid references with the highest correlation values are returned.
#' @details The map plotting of this function makes use of the functions of the \code{maps} package.
#' This method also makes use of the \code{\link[stats]{cor.test}} function from the \code{stats} package. When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
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
#' @author Ardern Hulme-Beaman
#'
#' @examples
#' Range.Exp <- .5
#'
#' Long.Range <- c(80, -100)
#' #Lat.Range <- c(floor(min(RatVarData$Lat))-Range.Exp,ceiling(max(RatVarData$Lat)+Range.Exp))
#'
#' #R.Samp <- c(30, 50)
#'
#'
#' #VarianceOrigin(LatLongs = cbind(RatVarData$Lat, RatVarData$Long),
#'  #              VarData = RatVarData$mean.of.sq.prodist.var,
#'   #             LongRange = Long.Range,
#'    #            LatRange = Lat.Range,
#'     #           RangeSamp = R.Samp,
#'      #          Method = 'Pearson',
#'       #         PacificCent = TRUE)
#'
#' #points(x = Rpraetor$Lat.Long$Long[1], y=Rpraetor$Lat.Long$Lat[1], col='blue', pch=16)
#'
#' @export
#'
#' @import mapdata


HeatMapping <- function(LatLongs, VarData, LongRange=c(0,0), LatRange=c(0,0), RangeSamp=10, Verbose=TRUE, PrintProg=FALSE, PlotRes=TRUE, HeatHue= c(.15, 1), TileSize=2, Method=c('Spearman', 'Pearson'), PacificCent=FALSE){

  #LatLongs=cbind(RatVarData$Lat, RatVarData$Long)

  #VarData=RatVarData$mean.of.sq.prodist.var
  #Range.Exp <- .5

  #LongRange <- c(75, -100)
  #LatRange <- c(floor(min(RatVarData$Lat))-Range.Exp,ceiling(max(RatVarData$Lat)+Range.Exp))
  #RangeSamp <- c(12, 42)
  #PacificCent=TRUE



  UserInputAssessment(LatLongs, RefData=NULL, Method = Method, RefDistMat = NULL, DistVec = NULL)


  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")


  #creating an empty object to be populated by results
  CoordsHeat <- NULL


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
  if (PacificCent==TRUE){
    MidRange <- seq(LongRange[1], LongRange[2]+360, by = LongRangeSteps*-1)
    #MidRange[which(MidRange>=180)] <- MidRange[which(MidRange>=180)]-360
    Longways <- c(LongRange[1]+LongRangeSteps, MidRange, LongRange[2]+(LongRangeSteps*-1)+360)
  } else {
    Longways <- c(LongRange[1]-LongRangeSteps, seq(LongRange[1], LongRange[2], by = LongRangeSteps), LongRange[2]+LongRangeSteps)
  }

  Latways <- c(LatRange[1]-LatRangeSteps, seq(LatRange[1], LatRange[2], by = LatRangeSteps), LatRange[2]+LatRangeSteps)


  #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
  for (i in Longways){
    for (j in Latways){
      #i <- Longways[1]
      #j <- Latways[1]
      coord <- c(j, i)

      GeographicDist <- suppressWarnings(GeoDist2Point(RefLatLongs = LatLongs[,], TargetLatLong = coord))

      if (Method=='Spearman'){
        #running the correlation to generate r
        CorRes <- suppressWarnings(stats::cor.test(x = VarData, y = GeographicDist, method = "spearman"))
      } else if (Method=='Pearson'){
        CorRes <- suppressWarnings(stats::cor.test(x = VarData, y = GeographicDist, method = "pearson"))
      }

      #combining the individual results and organising them
      results <- c(as.character(j),i, CorRes$estimate)

      #adding the results of this round to the previous results
      CoordsHeat <- rbind(CoordsHeat, results)
    }

    if (PrintProg==TRUE){
      svMisc::progress(value = which(Longways==i), max.value = length(Longways), progress.bar = FALSE)
      Sys.sleep(0.01)

    }
  }



  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat, row.names = 1:dim(CoordsHeat)[1])

  #naming the variables
  names(CoordsHeat) <- c("Lats", "Longs", "Cor")

  if (PacificCent==TRUE){
    PlottingMap <- 'mapdata::world2Hires'
  } else {
    PlottingMap <- "world"
  }

  if (PlotRes==TRUE){
    maps::map(PlottingMap, xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"))


    #creating colour scale from max and min correlation based on which variable we're using
    CoordsHeatNum <- chr2nu(CoordsHeat$Cor)
    OriginLoc <- CoordsHeat[which(CoordsHeatNum==max(CoordsHeatNum)),]


    CoordsHeatscaled <- (((CoordsHeatNum-min(CoordsHeatNum))/(max(CoordsHeatNum)-min(CoordsHeatNum)))*-1)+1

    CoordsHeats <- grDevices::hsv(h = HeatHue[1], v = 1, s = CoordsHeatscaled, alpha = HeatHue[2])

    #plotting the correlations
    graphics::points(x = as.character(CoordsHeat$Longs), y = as.character(CoordsHeat$Lat), pch=15, col=CoordsHeats,  cex=TileSize)


    maps::map(PlottingMap, xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"), add=T)

    graphics::points(x = as.character(LatLongs$Longs), y = as.character(LatLongs$Lat), pch=23, bg='orange',  cex=1)
    maps::map.axes()
  }

  if (Verbose==TRUE){
    return(CoordsHeat)
  } else {
    OriginLocCor <- CoordsHeat[which(CoordsHeat$Cor==max(chr2nu(CoordsHeat$Cor), na.rm = TRUE)),]

    return(list(Cor=OriginLocCor))
  }


}







#' Heat mapping function to identify region of highest value with significance
#'
#' This function takes a metric for spatial specimens or populations with defined locations
#' and calculates the region of highest values. This can be used with population variance
#' scores to identify regions of highest variance. A permutation test accompanies this to
#' identify regions where the correlation between the metric and geographic distance is
#' significant to the user determined level. This function runs permutations in parallel.
#' @param Perm an integer of the number of permutations the function is to run.
#' @param SigLevel a numeric value between 1 and 0 to set the required significance level for plotting a significant ploygon. Default is set to 0.05.
#' @inheritParams HeatMapping
#' @return If Verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If Verbose is FALSE then only those grid references with the highest correlation values are returned.
#' @details The map plotting of this function makes use of the functions of the \code{maps} package.
#' This method also makes use of the \code{\link[stats]{cor.test}} function from the \code{stats} package. When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
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
#' @author Ardern Hulme-Beaman
#'
#' @examples
#' Range.Exp <- .5
#'
#' Long.Range <- c(80, -100)
#' #Lat.Range <- c(floor(min(RatVarData$Lat))-Range.Exp,ceiling(max(RatVarData$Lat)+Range.Exp))
#'
#' #R.Samp <- c(30, 50)
#'
#'
#' #VarianceOrigin(LatLongs = cbind(RatVarData$Lat, RatVarData$Long),
#'  #              VarData = RatVarData$mean.of.sq.prodist.var,
#'   #             LongRange = Long.Range,
#'    #            LatRange = Lat.Range,
#'     #           RangeSamp = R.Samp,
#'      #          Method = 'Pearson',
#'       #         PacificCent = TRUE)
#'
#' #points(x = Rpraetor$Lat.Long$Long[1], y=Rpraetor$Lat.Long$Lat[1], col='blue', pch=16)
#'
#' @export


#RatVarData <- read.csv('C:/Users/Arder/Desktop/R work/Package Building/DetailedReg_Var_inter_photographer_Re_for_AR.csv')
#VarianceOrigin(LatLongs = cbind(RatVarData$Lat, RatVarData$Long),VarData = RatVarData$mean.of.sq.prodist.var,LongRange = Long.Range,LatRange = Lat.Range,RangeSamp = R.Samp,Method = 'Pearson',PacificCent = TRUE)
#LatLongs = cbind(RatVarData$Lat, RatVarData$Long);VarData = RatVarData$mean.of.sq.prodist.var;LongRange = Long.Range;LatRange = Lat.Range;RangeSamp = R.Samp;Method = 'Pearson';PacificCent = TRUE; Perm=10

HeatMappingPerm <- function(LatLongs, VarData, LongRange=c(0,0), LatRange=c(0,0), RangeSamp=10, Verbose=TRUE, PlotRes=TRUE, HeatHue= c(.15, 1), TileSize=2, Method=c('Spearman', 'Pearson'), PacificCent=FALSE, Perm=100, SigLevel=0.05){

  UserInputAssessment(LatLongs, RefData=NULL, Method, RefDistMat = NULL, DistVec = NULL)



  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")


  #creating an empty object to be populated by results
  CoordsHeat <- NULL


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
  if (PacificCent==TRUE){
    MidRange <- seq(LongRange[1], LongRange[2]+360, by = LongRangeSteps*-1)
    #MidRange[which(MidRange>=180)] <- MidRange[which(MidRange>=180)]-360
    Longways <- c(LongRange[1]+LongRangeSteps, MidRange, LongRange[2]+(LongRangeSteps*-1)+360)
  } else {
    Longways <- c(LongRange[1]-LongRangeSteps, seq(LongRange[1], LongRange[2], by = LongRangeSteps), LongRange[2]+LongRangeSteps)
  }

  Latways <- c(LatRange[1]-LatRangeSteps, seq(LatRange[1], LatRange[2], by = LatRangeSteps), LatRange[2]+LatRangeSteps)


  GeoDist2PointPar <- function(RefLatLongs, TargetLatLong){
    GeographicDist <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1], r=6378.137)
    return(GeographicDist)
  }

  VarianceOriginPar <- function(LatLongs, VarData, Longways, Latways, Verbose, Method, PacificCent){
    VarDataPerm <- VarData[sample(1:length(VarData), size = length(VarData))]
    CoordsHeatPar <- matrix(NA, nrow = length(Latways), ncol = length(Longways))

    #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
    for (i in 1:length(Longways)){
      for (j in 1:length(Latways)){
        #i <- j <- 1
        I <- Longways[i]
        J <- Latways[j]
        coord <- c(J, I)

        GeographicDist <- suppressWarnings(GeoDist2PointPar(RefLatLongs = LatLongs[,], TargetLatLong = coord))

        if (Method=='Spearman'){
          #running the correlation to generate r
          CorRes <- suppressWarnings(stats::cor.test(x = VarDataPerm, y = GeographicDist, method = "spearman"))
        } else if (Method=='Pearson'){
          CorRes <- suppressWarnings(stats::cor.test(x = VarDataPerm, y = GeographicDist, method = "pearson"))
        }

        CoordsHeatPar[j,i] <- CorRes$estimate

      }

    }


    return(CoordsHeatPar)

  }

  arraybind <- function(A, a){
    #A=PermVarHeat; a=testRes

    if (is.null(A)){
      return(a)
    }

    array2mat <- function(Array){

      Matrix <- matrix(NA, nrow = dim(Array)[3], ncol = length(c(t(Array[,,1]))))
      for (i in 1:dim(Array)[3]){
        #i <- 1
        Matrix[i,] <- c(t(Array[,,i]))
      }
      return(Matrix)
    }

    mat2array <- function(mat, arraydim){
      NewArray <- array(data = NA, dim = c(dim(mat)[2]/arraydim, arraydim, dim(mat)[1]))


      for (i in 1:dim(mat)[1]){
        #i <- 1
        Mat4Array <- matrix(as.numeric(mat[i,]), nrow = dim(mat)[2]/arraydim, ncol = arraydim, byrow = TRUE)
        NewArray[,,i] <- Mat4Array
      }
      return(NewArray)

    }

    if (is.matrix(A)){
      if (is.matrix(a)){
        combinedData <- rbind(c(A), c(a))
      }
    } else if (is.matrix(a)){
      combinedData <- rbind(array2mat(A), c(a))
    } else {
      combinedData <- rbind(array2mat(A), array2mat(a))
    }


    return(mat2array(combinedData, arraydim = dim(A)[2]))


  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)


  #Perm=1000
  a <- 1
  PermVarHeat <- foreach::foreach(a = 1:Perm, .combine = arraybind, .packages='GeoOrigins') %dopar% {
    VarianceOriginPar(LatLongs, VarData, Longways, Latways, Verbose, Method, PacificCent)
  }

  parallel::stopCluster(clust)

  TrueVarVals <- HeatMapping(LatLongs, VarData, LongRange, LatRange, RangeSamp, Verbose=TRUE, PlotRes=FALSE, HeatHue=HeatHue, TileSize=TileSize, Method = Method, PacificCent, PrintProg = FALSE)

  if (PacificCent==TRUE){
    PlottingMap <- 'mapdata::world2Hires'
    Longways[which(Longways<=0)] <- Longways[which(Longways<=0)]+360
  } else {
    PlottingMap <- "world"
  }


  maps::map(PlottingMap, xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"))


  #creating colour scale from max and min correlation based on which variable we're using
  TrueVarValsNum <- chr2nu(TrueVarVals$Cor)

  TrueVarValsscaled <- (((TrueVarValsNum-min(TrueVarValsNum))/(max(TrueVarValsNum)-min(TrueVarValsNum)))*-1)+1

  TrueVarValss <- grDevices::hsv(h = HeatHue[1], v = 1, s = TrueVarValsscaled, alpha = HeatHue[2])

  #plotting the correlations
  graphics::points(x = as.character(TrueVarVals$Longs), y = as.character(TrueVarVals$Lat), pch=15, col=TrueVarValss,  cex=TileSize)

  maps::map(PlottingMap, xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"), add = TRUE)

  graphics::points(x = as.character(LatLongs$Longs), y = as.character(LatLongs$Lat), pch=23, bg='orange',  cex=1)
  maps::map.axes()

  TrueVarValsMat <- matrix(as.numeric(as.character(TrueVarVals$Cor)), nrow = length(Latways), ncol = length(Longways))

  SigCells <- PermArraySig(PermArray=PermVarHeat, TrueMatrix=TrueVarValsMat, percentile = SigLevel)
  SigPlot <- TrueVarVals[which(SigCells==TRUE),]

  contour.95 <-  Construct_contour(SigPlot[,1:2])
  polycol <- grDevices::hsv(h = HeatHue[1], s = 1, v = .8, alpha = HeatHue[2])

  graphics::polygon(contour.95$x, contour.95$y, col=NA, border=polycol, lwd=2)

  dimnames(PermVarHeat) <- list(Latways, Longways, 1:Perm)

  VarResults <- list(TrueVarCorValues=TrueVarVals, PermutationArray=PermVarHeat)

  return(VarResults)
}




#' Two location permutation significance signficance test
#'
#' This function takes the output from the \code{HeatMappingPerm} function and
#' two user selected locations and assesses whether there is a significant difference
#' between the identified correlation values at the two locations. The method does
#' not look specifically at the user defined locations, but rather looks at the nearest
#' grid square in the \code{HeatMappingPerm} output; this is primarily for computing
#' efficiency.
#' @param LocationA a numeric vector detailing the Latitude and Longitude values of the first location the user selects.
#' @param LocationB a numeric vector detailing the Latitude and Longitude values of the second location the user selects.
#' @param HeatMappingPermResults the output from the \code{HeatMappingPerm} function. This should be a list of two objects: TrueVarCorValues and PermutationArray.
#' @param SigLevel a numeric value between 1 and 0 to set the required significance level for plotting a significant ploygon. Default is set to 0.05.
#' @inheritParams HeatMapping
#' @return The difference between the two locations and what level of significance this is.
#'
#' @author Ardern Hulme-Beaman
#'


TwoLocationTest <- function(LocationA, LocationB, HeatMappingPermResults){
  #LocationA=c(0, 90)
  #LocationB=c(0, 140)
  #HeatMappingPermResults=VarResults


  GridLocations <- HeatMappingPermResults$TrueVarCorValues[,-3]

  Dist2A <- suppressWarnings(GeoDist2Point(cbind(chr2nu(GridLocations$Lats), chr2nu(GridLocations$Longs)), LocationA))
  SelectedLocationAPos <- which(Dist2A==min(Dist2A))
  SelectedLocationA <- GridLocations[SelectedLocationAPos,]
  SelectedLocationAVal <- HeatMappingPermResults$TrueVarCorValues[SelectedLocationAPos,3]

  returnedLocA <- rbind(LocationA, SelectedLocationA)
  rownames(returnedLocA) <- c('UserSelectedLocation', 'NearestGridLocation')


  ALongIndex <- which(colnames(HeatMappingPermResults$PermutationArray)==as.character(SelectedLocationA[2]))
  ALatIndex <- which(rownames(HeatMappingPermResults$PermutationArray)==as.character(SelectedLocationA[1]))


  Dist2B <- suppressWarnings(GeoDist2Point(cbind(chr2nu(GridLocations$Lats), chr2nu(GridLocations$Longs)), LocationB))
  SelectedLocationBPos <- which(Dist2B==min(Dist2B))
  SelectedLocationB <- GridLocations[SelectedLocationBPos,]
  SelectedLocationBVal <- HeatMappingPermResults$TrueVarCorValues[SelectedLocationBPos,3]

  returnedLocB <- rbind(LocationB, SelectedLocationB)
  rownames(returnedLocB) <- c('UserSelectedLocation', 'NearestGridLocation')


  BLongIndex <- which(colnames(HeatMappingPermResults$PermutationArray)==as.character(SelectedLocationB[2]))
  BLatIndex <- which(rownames(HeatMappingPermResults$PermutationArray)==as.character(SelectedLocationB[1]))


  TrueDist <- dist(c(SelectedLocationBVal, SelectedLocationAVal))

  ALocVals <- HeatMappingPermResults$PermutationArray[ALatIndex, ALongIndex,]
  BLocVals <- HeatMappingPermResults$PermutationArray[BLatIndex, BLongIndex,]

  PerDists <- apply(rbind(ALocVals, BLocVals), MARGIN = 2, FUN = dist)

  #hist(PerDists)
  Sig <- sum(PerDists>TrueDist)/length(PerDists)

  CorVals <- as.numeric(c(SelectedLocationAVal, SelectedLocationBVal))
  Locations <- c('Location A', 'Location B')
  minVal <- which(CorVals==min(CorVals))


  PrintOut <- paste(Locations[minVal], ' ', round(as.numeric(SelectedLocationAVal), 4), ' is signficantly less than ', Locations[-minVal],' ', round(as.numeric(SelectedLocationBVal), 4), ' with an approximate p value of ', Sig, sep='')
  print(PrintOut)

  return(list(LocationA = returnedLocA, LocationB = returnedLocB, LocationAValue=as.numeric(SelectedLocationAVal), LocationBValue=as.numeric(SelectedLocationBVal), PermutatedDistribution=PerDists))

}

