

#' Wrapper for procdist function to output a distance table
#'
#' This function builds a square matrix of pairwise procrustes distances among specimens using
#' the \code{procdist} function from the \code{shapes} package.
#' @param A is a three dimensional data array of landmark x landmark dimensions x specimens
#' @return This function returns a square matrix of Procrustes distances, which is required for both the \code{IDbyDistance.RawData.CCV} and the \code{BoundaryFinder} functions.
#' @section Citations:
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#'
#' @export



ProcDistanceTable <- function(A){
  #A <- GPA$coords

  ProcDtable <- matrix(0, dim(A)[3], dim(A)[3])
  for (i in 1:dim(A)[3]){
    for (j in 1:dim(A)[3]){
      if (j<i){
        #i <- 1
        #j <- 2
        ProcDtable[i,j] <- shapes::procdist(A[,,i], A[,,j], type="full")
      }
    }
  }

  ProcDrableRes <- as.matrix(as.dist(ProcDtable))
  return(ProcDrableRes)
}


