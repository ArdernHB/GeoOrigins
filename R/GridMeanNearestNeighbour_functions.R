


#' Spatially provenance a specimen to a grid cell by nearest neighbour methods
#'
#' This function takes the raw variables of an unknown specimen and reference specimens
#' and uses euclidean distances to calculate a likely spatial provenance. Note that this
#' procedure can only be applied to one unknown specimen at a time. Shape variables
#' can be specified and if so Procrustes distances can be calculated. It has two applications
#' either: calculating a specimens' provenance, or alternatively it can be used to calculate the
#' minimum correlation coefficient needed to correctly identify a known specimen at its true
#' collection location. The second application of this function can work as a correct cross-
#' validation process if applied in a for loop.
#' @inheritParams IDbyDistanceRawData
#' @inheritParams BoundaryFinder
#'
#' @details When used for shape data and for Procrustes distances this function makes use of the \code{\link[shapes]{procGPA}} and \code{\link[shapes]{procdist}} functions from the \code{shapes} package. When Euclidean distances are employed the \code{\link[stats]{dist}} function of the base \code{stats} package is used.
#' This method also makes use of the \code{\link[stats]{cor.test}} function from the \code{stats} package. When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
#' The map plotting of this function makes use of the functions of the \code{maps} package.
#'
#' @section Citations:
#'
#' Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
#' Enhancements by Thomas P Minka and Alex Deckmyn. (2017). maps: Draw Geographical Maps. R
#' package version 3.2.0. https://CRAN.R-project.org/package=maps
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#' Grosjean, Ph. (2016). svMisc: SciViews-R. UMONS, Mons, Belgium.
#' http://www.sciviews.org/SciViews-R.
#'
#' @author Ardern Hulme-Beaman
#'
#' @examples
#' Range.Exp <- .5
#'
#' Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Long)+Range.Exp))
#' Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Lat)+Range.Exp))
#'
#' RpraetorDataMat <- Array2Mat(Rpraetor$LMs)
#' R.Samp <- c(7,15)
#'
#' GridIDRes <- NULL
#'
#' #for (i in 1:dim(RpraetorDataMat)[1]){
#'
#'
#' #iRes <- IDbyGridMeanRawData(LatLongs = Rpraetor$Lat.Long[-i,],
#'  #                        TargetData = RpraetorDataMat[i,],
#'  #                        RefData = RpraetorDataMat[-i,],
#'  #                        ShapeData = TRUE,
#'  #                        ShapeDim = 2,
#'  #                        LongRange = Long.Range,
#'  #                        LatRange = Lat.Range,
#'  #                        RangeSamp = R.Samp,
#'  #                        Validate = TRUE,
#'  #                        ValidLatLongs = Rpraetor$Lat.Long[i,],
#'  #                        IgnorePrompts = TRUE,
#'  #                        PlotProv = FALSE)
#'
#'   #GridIDRes <- rbind(GridIDRes, c(as.matrix(iRes)[1,], dim(iRes)[1]))
#' #}
#'
#' #table(GridIDRes[,5])/length(GridIDRes[,5])
#'

IDbyGridMeanRawData <- function(LatLongs, TargetData, RefData, ShapeData=TRUE, ShapeDim=2, LongRange=c(0,0), LatRange=c(0,0), ExpandMap=c(0,0), RangeSamp=10, PrintProg=FALSE, Validate=FALSE, ValidLatLongs, HeatHue= c(.15, 1), TileSize=2, IgnorePrompts=FALSE, PlotProv=TRUE){

  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")



  #organising data for ease of analysis
  #combining ref and target specimens with target first
  if (ShapeData==TRUE){
    TotalShape.raw <- rbind(TargetData, RefData)
    gpaRes <- shapes::procGPA(Mat2Array(TotalShape.raw, LMdim = ShapeDim))
    TotalShape <- Array2Mat(gpaRes$rotated)
  } else {
    TotalShape <- rbind(TargetData, RefData)
  }


  #creating a range that will cover the whole geographic area of interest
  #this is for the function to loop through later on

  UserGood <- 'bad'
  while (UserGood!='good') {

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

    LongSpace <- c(Longways[2]-Longways[1])/2
    LatSpace <- c(Latways[2]-Latways[1])/2

    GridLinesLong <- Longways+LongSpace
    GridLinesLat <- Latways+LatSpace

    if (PlotProv==TRUE){
      maps::map("world", xlim=c(min(Longways)-ExpandMap[2], max(Longways)+ExpandMap[2]), ylim=c(min(Latways)-ExpandMap[1], max(Latways)+ExpandMap[1]), interior=FALSE, col="black", bg=graphics::par(bg="white"))
      graphics::points(x = as.character(LatLongs$Long), y = as.character(LatLongs$Lat), pch=23, bg='blue',  cex=1)
      graphics::abline(h=GridLinesLat, v = GridLinesLong)
      maps::map.axes()

      #GridCellCentres <- expand.grid(Longways, Latways)
      #points(GridCellCentres)

      if (IgnorePrompts==FALSE){

        UserInput <- readline(prompt = "Do you wish to continue with this grid division? Please respond y or n")

        while (!(UserInput%in% c('n','y'))){
          UserInput <- readline(prompt = "Response was neither n or y. Do you wish to continue with this grid division? Please respond y or n")
        }

        if (UserInput=='n'){
          RangeSampLong <- readline('What number of longitudinal divisions do you want?')
          RangeSampLat <- readline('What number of latitudinal divisions do you want?')
          if (is.na(suppressWarnings(as.numeric(RangeSampLong))) || is.na(suppressWarnings(as.numeric(RangeSampLat)))){
            stop('Error: user provided a non numeric value. Please rerun function and adjust the value of the RangeSamp argument')
          } else {
            RangeSamp <- as.numeric(c(RangeSampLat, RangeSampLong))

          }
        } else if (UserInput=='y'){
          UserGood <- 'good'
        }

      } else {UserGood <- 'good'}

    } else {
      UserGood <- 'good'
    }



  }


  #creating an empty object to be populated by results
  CoordsHeat <- NULL


  #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
  for (i in 1:length(Longways)){
    for (j in 1:length(Latways)){
      #i <- 2
      #j <- 6
      coord <- c(Latways[j], Longways[i])

      GridCellSelect <- which(LatLongs$Lats>c(coord[1]-LatSpace) & LatLongs$Lats<c(coord[1]+LatSpace) & LatLongs$Longs<c(coord[2]+LongSpace) &  LatLongs$Longs>c(coord[2]-LongSpace))
      #points(x = LatLongs[GridCellSelect,2], y = LatLongs[GridCellSelect,1], pch=18, col='red')

      if (length(GridCellSelect)>0){
        if (ShapeData==TRUE && length(GridCellSelect)>1){
          CellGPA <- shapes::procGPA(Mat2Array(RefData[GridCellSelect,], ShapeDim))

          Dist2CellMean <- shapes::procdist(CellGPA$mshape, matrix(TargetData, nrow = length(TargetData)/ShapeDim, ncol = ShapeDim, byrow = TRUE), type = 'full')
        } else if (ShapeData==TRUE && length(GridCellSelect)==1){
          Dist2CellMean <- shapes::procdist(matrix(RefData[GridCellSelect,], nrow = length(TargetData)/ShapeDim, ncol = ShapeDim, byrow = TRUE),
                                            matrix(TargetData, nrow = length(TargetData)/ShapeDim, ncol = ShapeDim, byrow = TRUE), type = 'full')
        } else if (ShapeData==FALSE && length(GridCellSelect)>1){
          Dist2CellMean <- stats::dist(rbind(TargetData, colMeans(RefData[GridCellSelect,])))

        } else if (ShapeData==FALSE && length(GridCellSelect)==1){
          Dist2CellMean <- stats::dist(rbind(TargetData, RefData[GridCellSelect,]))

        }

        CoordsHeat <- rbind(CoordsHeat, c(coord, Dist2CellMean))


      }
      if (PrintProg==TRUE){
        svMisc::progress(value = which(Longways==Longways[i]), max.value = length(Longways), progress.bar = FALSE)
        #Sys.sleep(0.01)

      }
    }
  }

  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat)

  #naming the variables
  names(CoordsHeat) <- c("Lat", "Long", "Trait.Dist.to.Unknown")

  CoordsHeatNum <- chr2nu(CoordsHeat$Trait.Dist.to.Unknown)
  OriginLoc <- CoordsHeat[which(CoordsHeatNum==min(CoordsHeatNum)),]

  if (PlotProv==TRUE){
    CoordsHeatscaled <- 1-(CoordsHeatNum-min(CoordsHeatNum))/(max(CoordsHeatNum)-min(CoordsHeatNum))

    CoordsHeats <- grDevices::hsv(h = HeatHue[1], v = 1, s = CoordsHeatscaled, alpha = HeatHue[2])

    graphics::points(x = as.character(CoordsHeat$Long), y = as.character(CoordsHeat$Lat), pch=22, bg=CoordsHeats,  cex=TileSize)
    graphics::points(x = as.character(OriginLoc$Long), y = as.character(OriginLoc$Lat), pch=0, col='red',  cex=TileSize, lwd=2)

  }

  if (Validate==TRUE){
    #ValidLatLongs <- Rpraetor$Lat.Long[1,]
    TrueGridComb <- rbind(c(as.numeric(ValidLatLongs), 0), CoordsHeat)
    Dist2True <- GeoDist2Point(RefLatLongs = TrueGridComb[-1,], TargetLatLong = TrueGridComb[1,])

    SortedDist2True <- sort(Dist2True, index.return=TRUE)


    SortedGridMeanDistResPrep <- cbind(CoordsHeat[SortedDist2True$ix,], SortedDist2True$x)

    if (dim(OriginLoc)[1]>1){
      SortedGridMeanDistRes <- cbind(SortedGridMeanDistResPrep, c('incorrect', rep('', length(SortedGridMeanDistResPrep[,1])-1)))
    } else if (sum(SortedGridMeanDistResPrep[1,1:3]==OriginLoc)==3){
      SortedGridMeanDistRes <- cbind(SortedGridMeanDistResPrep, c('correct', rep('', length(SortedGridMeanDistResPrep[,1])-1)))
    } else {
      SortedGridMeanDistRes <- cbind(SortedGridMeanDistResPrep, c('incorrect', rep('', length(SortedGridMeanDistResPrep[,1])-1)))
    }

    names(SortedGridMeanDistRes)[4:5] <- c('Geo.Dist.to.True.Location', "CCV")

  } else {

    SortedGridMeanDistRes <- CoordsHeat[sort(CoordsHeat$Trait.Dist.to.Unknown, index.return=TRUE)$ix,]
  }

  return(SortedGridMeanDistRes)

}


