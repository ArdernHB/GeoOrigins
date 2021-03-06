% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GridMeanNearestNeighbour_functions.R
\name{IDbyGridMeanRawData}
\alias{IDbyGridMeanRawData}
\title{Spatially provenance a specimen to a grid cell by nearest neighbour methods}
\usage{
IDbyGridMeanRawData(
  LatLongs,
  TargetData,
  RefData,
  ShapeData = TRUE,
  ShapeDim = 2,
  LongRange = c(0, 0),
  LatRange = c(0, 0),
  ExpandMap = c(0, 0),
  RangeSamp = 10,
  PrintProg = FALSE,
  Validate = FALSE,
  ValidLatLongs,
  HeatHue = c(0.15, 1),
  TileSize = 2,
  IgnorePrompts = FALSE,
  PlotProv = TRUE
)
}
\arguments{
\item{LatLongs}{a matrix of n rows by 2 columns where n is the number of reference specimens in your dataset and the columns are Latitude and Longitude values in that order. These latitude-longitude coordinates should be of the locations of the reference specimens.}

\item{TargetData}{is a vector of unknown specimen data. If it is geometric morphometric data is should be a vector of superimposed coordinated in the format X1, Y1, X2, Y2... etc. (or a vector of other standardised variables that can appropriately have euclidean distances calculated between it and reference variables). NB if applied to a reference collection specimen be sure to remove it from the RefData dataset.}

\item{RefData}{is a matrix of Reference specimen data where the rows are the individual reference specimens and the columns are the variables in the same order as the TargetData vector.}

\item{ShapeData}{logical indicating whether the data is geometric morphometric shape data that requires superimposition. Default set to TRUE}

\item{ShapeDim}{integer either 2 or 3 to indicate the dimensions of landmark coordinates if the data is geometric morphometric data.}

\item{LongRange}{is a vector of 2 elements defining the maximum and minimum Longitude values that the provenancing method should explore. This will also define the mapping range in the final plotted output.}

\item{LatRange}{is a vector of 2 elements defining the maximum and minimum Latitude values that the provenancing method should explore. This will also define the mapping range in the final plotted output.}

\item{ExpandMap}{is a vector of 2 elements for expanding the plotting region of the map. The first element expands the latitudinal area and the second element expands the longitudinal area.}

\item{RangeSamp}{is an integer vector of 1 or 2 elements that defines the resolution of spatial sampling. If one element is provided then both the latitude and longitude ranges are equally and evenly sampled using this value. If 2 elements are provided they should be in the order of latitude longitude and each range will be evenly sampled with its respective value.}

\item{PrintProg}{logical whether or not to print a progress bar. Default set to FALSE.}

\item{Validate}{logical whether or not to run a correct cross-validation analysis to find the lowest required correlation value for correct identification.}

\item{ValidLatLongs}{if the process is carried out on a specimen of known location `(e.g. Validate=TRUE)`, then the latitude longitude coordinates for that location should be provided here in that order.}

\item{HeatHue}{numeric vector of 2 elements each between 0 and 1. The first should be the hue value on the HSV scale; the second value should be the level of transparency of the colour used.}

\item{TileSize}{numeric to dictate the pixel size of the heat mapping colour values.}

\item{IgnorePrompts}{default is set to FALSE, but if set to TRUE queries such as those confirming the location of the datadump will be suppressed.}

\item{PlotProv}{logical if the map should be printed with a polygon demarking a contour at a user defined correlation value.}
}
\description{
This function was for comparison between a geographically spatially based k-NN approach
compared with the GeoOrigins correlation approach; it is the closest direct comparison between
k-NN and the new correlation method. However, it is still effectively a user defined group based method
because the user defines the grid size and the members of each grid cell effectively form a group
that unknowns are identified to (albeit that the unknowns are identified not to the group as such but the
group mean i.e. the mean of all the specimens in the grid cell).
This function takes the raw variables of an unknown specimen and reference specimens
and uses euclidean distances to calculate a likely spatial provenance. Note that this
procedure can only be applied to one unknown specimen at a time. Shape variables
can be specified and if so Procrustes distances can be calculated. It has two applications
either: calculating a specimens' provenance, or alternatively it can be used to calculate the
minimum correlation coefficient needed to correctly identify a known specimen at its true
collection location. The second application of this function can work as a correct cross-
validation process if applied in a for loop.
}
\details{
When used for shape data and for Procrustes distances this function makes use of the \code{\link[shapes]{procGPA}} and \code{\link[shapes]{procdist}} functions from the \code{shapes} package. When Euclidean distances are employed the \code{\link[stats]{dist}} function of the base \code{stats} package is used.
This method also makes use of the \code{\link[stats]{cor.test}} function from the \code{stats} package. When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
The map plotting of this function makes use of the functions of the \code{maps} package.
}
\section{Citations}{


Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
Enhancements by Thomas P Minka and Alex Deckmyn. (2017). maps: Draw Geographical Maps. R
package version 3.2.0. https://CRAN.R-project.org/package=maps

Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
https://CRAN.R-project.org/package=shapes

Grosjean, Ph. (2016). svMisc: SciViews-R. UMONS, Mons, Belgium.
http://www.sciviews.org/SciViews-R.
}

\examples{
Range.Exp <- .5

Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long))
              -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Long)+Range.Exp))
Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat))
              -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Lat)+Range.Exp))

RpraetorDataMat <- Array2Mat(Rpraetor$LMs)
R.Samp <- c(7,15)

GridIDRes <- NULL

#for (i in 1:dim(RpraetorDataMat)[1]){


#iRes <- IDbyGridMeanRawData(LatLongs = Rpraetor$Lat.Long[-i,],
 #                        TargetData = RpraetorDataMat[i,],
 #                        RefData = RpraetorDataMat[-i,],
 #                        ShapeData = TRUE,
 #                        ShapeDim = 2,
 #                        LongRange = Long.Range,
 #                        LatRange = Lat.Range,
 #                        RangeSamp = R.Samp,
 #                        Validate = TRUE,
 #                        ValidLatLongs = Rpraetor$Lat.Long[i,],
 #                        IgnorePrompts = TRUE,
 #                        PlotProv = FALSE)

  #GridIDRes <- rbind(GridIDRes, c(as.matrix(iRes)[1,], dim(iRes)[1]))
#}

#table(GridIDRes[,5])/length(GridIDRes[,5])

}
\author{
Ardern Hulme-Beaman
}
\keyword{internal}
