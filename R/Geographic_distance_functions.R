

#' Returns vector of Haversine distances
#'
#' This function uses the \code{geosphere} package function \code{\link[geosphere]{distHaversine}} to calculate multiple distances between a set of Latitude Longitude coordinates to a target Latitude Longitude coordinate.
#' @param RefLatLongs a matrix of n rows by 2 columns where n is the number of coordinates and the columns are Latitude and Longitude values in that order.
#' @param TargetLatLong a vector of 2 elements that are Latitude and Longitude values in that order. This is the coordinate of interest which all distances will be calculated to.
#' @return The distance in metres between the point of interest (TargetLatLong) to all other points in the reference material (RefLatLongs).
#' @section Citations:
#'
#' Sinnott, R.W, 1984. Virtues of the Haversine. Sky and Telescope 68(2): 159
#'
#' Robert J. Hijmans (2017). geosphere: Spherical Trigonometry. R package version 1.5-7.
#' https://CRAN.R-project.org/package=geosphere
#'
#' @keywords internal
#' @keywords geographic distance
#' @keywords Haversine distance
#' @author Ardern Hulme-Beaman
#' @import geosphere


GeoDist2Point <- function(RefLatLongs, TargetLatLong){

  GeographicDist <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1], r=6378.137)

  return(GeographicDist)
}



#' Returns a table of Haversine distances
#'
#' This function uses the \code{geosphere} package function \code{\link[geosphere]{distHaversine}} to calculate pairwise distances among all inputted Latitude Longitude coordinates.
#' @param RefLatLongs a matrix of n rows by 2 columns where n is the number of coordinates and the columns are Latitude and Longitude values in that order.
#' @param IDs a vector of unique IDs that correspond with (and are in the same order as) the Latitude Longitude coordinates. These IDs will then be used to name the columns and the rows of the returned table of pairwise distances. The default is set to NA. If set to NA the returned table will not have named columns or rows.
#' @return A square matrix of pairwise distances in metres among all inputted coordinates (RefLatLongs).
#' @section Citations:
#'
#' Sinnott, R.W, 1984. Virtues of the Haversine. Sky and Telescope 68(2): 159
#'
#' Robert J. Hijmans (2017). geosphere: Spherical Trigonometry. R package version 1.5-7.
#' https://CRAN.R-project.org/package=geosphere
#'
#' @keywords geographic distances
#' @keywords Haversine distances
#' @author Ardern Hulme-Beaman
#' @import geosphere


GeoDist2PointTable <- function(RefLatLongs, IDs=NA){

  CombinedLatLongs <- paste(RefLatLongs[,1], RefLatLongs[,2], sep = "_")
  CombinedLatLongsFacts <- as.factor(CombinedLatLongs)


  UniqueLatLongsCombined <- sort(unique(CombinedLatLongs))
  UniqueLatLongs <- do.call('rbind', strsplit(UniqueLatLongsCombined, split = "_"))

  GeographicDist <- matrix(NA, nrow = length(UniqueLatLongsCombined), ncol = length(UniqueLatLongsCombined))
  for (a in 1:length(UniqueLatLongsCombined)){
    for (b in 1:length(UniqueLatLongsCombined)){
      if (a>=b){
        distA <- geosphere::distHaversine(p1 = chr2nu(UniqueLatLongs[a,2:1]), p2 = chr2nu(UniqueLatLongs[b,2:1]), r=6378.137)
        GeographicDist[a,b] <- distA
      }
    }
  }

  GeographicDist.Mat <- as.matrix(stats::as.dist(GeographicDist))

  CompleteRes <- GeographicDist.Mat[CombinedLatLongsFacts, CombinedLatLongsFacts]

  if (length(IDs)<dim(RefLatLongs)[1] || length(IDs)==1 && is.na(IDs)==TRUE){
    return(CompleteRes)
  } else {
    colnames(CompleteRes) <- rownames(CompleteRes) <- IDs
    return(CompleteRes)
  }

}



