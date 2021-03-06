

#' Wrapper for procdist function to output a distance table
#'
#' This function builds a square matrix of pairwise procrustes distances among specimens using
#' the \code{\link[shapes]{procdist}} function from the \code{shapes} package.
#' @param A is a three dimensional data array of landmark x landmark dimensions x specimens
#' @return This function returns a square matrix of Procrustes distances, which is required for both the \code{IDbyDistanceRawDataCCV} and the \code{BoundaryFinder} functions.
#' @section Citations:
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#' @seealso
#' \code{\link{ProcDistanceTablePar}}
#'
#' @examples
#' RatDistMat <- ProcDistanceTable(Rpraetor$LMs)
#'
#' @import shapes
#' @export



ProcDistanceTable <- function(A){

  StartTime <- Sys.time()

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

  EndTime <- Sys.time()

  EndTime-StartTime

  ProcDrableRes <- as.matrix(stats::as.dist(ProcDtable))
  return(ProcDrableRes)
}






#' Wrapper for procdist function to output a distance table using parallel processing
#'
#' This function builds a square matrix of pairwise procrustes distances among specimens using
#' the \code{\link[shapes]{procdist}} function from the \code{shapes} package. This has been set
#' up to run parallel cores using the \code{foreach}, \code{doParallel} and \code{parallel} packages.
#' @param A is a three dimensional data array of landmark x landmark dimensions x specimens
#' @return This function returns a square matrix of Procrustes distances, which is required for both the \code{\link{IDbyDistanceRawDataCCV}} and the \code{BoundaryFinder} functions.
#' @section Citations:
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#' @seealso
#' \code{\link{ProcDistanceTable}}
#' @examples
#' #See ProcDistanceTable() function example
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import shapes
#' @export



ProcDistanceTablePar <- function(A){


  Res <- matrix(0, nrow = dim(A)[3], ncol = dim(A)[3])
  Res[lower.tri(Res)] <- 1
  Pairedindex <- which(Res==1, arr.ind = TRUE)


  DistMatloop <- function(X, index){
    ProcDRes <- shapes::procdist(X[,,index[1]], X[,,index[2]], type="full")
    return(c(index, ProcDRes))
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)


  a <- 1
  ProcDResComp <- foreach::foreach(a = 1:dim(Pairedindex)[1], .combine = rbind) %dopar%{
    DistMatloop(X=A, index = Pairedindex[a,])
  }

  parallel::stopCluster(clust)

  Res[lower.tri(Res)] <- as.numeric(ProcDResComp[,3])
  ProcDrableRes <- as.matrix(stats::as.dist(Res))
  return(ProcDrableRes)
}


