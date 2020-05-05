
#' Internal function for quick class change of a vector to numeric
#'
#' This function takes a factor vector and changes it to a numeric.
#' This is primarily to avoid instances where a csv file has been opened
#' without column classes being specified, which can result in numeric
#' values being accidently read as factors, which in turn when treated
#' with as.numeric() will accidently be given integer values in the order
#' of the factor levels. By treating the values as characters first it avoids
#' this potential error.
#' @return
#'
#' This function returns a numeric vector
#'
#' @param x is a vector to be converted to a numeric vector
#'
#' @keywords internal
#' @author Ardern Hulme-Beaman

chr2nu <- function(X){as.numeric(as.character(X))}





#' Simple function for conversion of a matrix of landmarks to an array
#'
#' This function takes a matrix where rows are specimens and columns are
#' alternating landmark coordinate values and converts it to a 3 dimensional
#' array where rows are landmarks and columns correspond with the dimensions.
#' @return
#'
#' This function returns an array
#'
#' @param Mat is a matrix of landmark data
#' @param LMdim is the number of dimensions of the data, either 2 or 3
#'
#' @keywords data format
#' @author Ardern Hulme-Beaman
#' @export



Mat2Array <- function(Mat, LMdim){
  NewArray <- array(data = NA, dim = c(dim(Mat)[2]/LMdim, LMdim, dim(Mat)[1]))


  for (i in 1:dim(Mat)[1]){
    #i <- 1
    Mat4Array <- matrix(as.numeric(Mat[i,]), nrow = dim(Mat)[2]/LMdim, ncol = LMdim, byrow = TRUE)
    NewArray[,,i] <- Mat4Array
  }


  return(NewArray)
}



#' Simple function for conversion of an array to a matrix
#'
#' This function takes an array and converts it to a matrix where rows
#' specimens and columns are alternating landmark variables
#' (e.g. X1, X2, Y1, Y2...)
#' @return
#'
#' This function returns an array
#'
#' @param Array is an array of landmark data
#'
#' @keywords data format
#' @author Ardern Hulme-Beaman
#' @export


Array2Mat <- function(Array){
  Matrix <- matrix(NA, nrow = dim(Array)[3], ncol = length(c(t(Array[,,1]))))
  for (i in 1:dim(Array)[3]){
    #i <- 1
    Matrix[i,] <- c(t(Array[,,i]))
  }
  return(Matrix)
}



#' Internal function for assessment of user input
#'
#' This function assesses whether user inputs are as expected and it not
#' the function throws an error
#'
#' @inheritParams GeoDist2Point
#' @inheritParams IDbyDistanceRawData
#' @inheritParams BoundaryFinder
#'
#' @author Ardern Hulme-Beaman




UserInputAssessment <- function(LatLongs, RefDistMat=NA, TargetData=NA, LongRange, LatRange, RangeSamp, Method){

  if (!(is.na(RefDistMat))){
    if (dim(LatLongs)[1]!=dim(RefDistMat)[1]){
      stop('Error: the number of latitude and longitude coordinates you have provided do not match the number of specimens in the dataset.
     \n Please check these match and most importantly are in the same order and rerun the function')
    }
  }

  if (!(is.na(TargetData))){
    if (dim(LatLongs)[1]!=dim(TargetData)[3]){
      stop('Error: the number of latitude and longitude coordinates you have provided do not match the number of specimens in the dataset.
     \n Please check these match and most importantly are in the same order and rerun the function')
    }
  }



  if (dim(LatLongs)[2]!=2){
    stop('Error: the matrix of latitude and longitude coordinates you have provided does not contain 2 columns.
     \n Please check that you have supplied the correct data and in the correct format.
     \n Importantly, please check that the supplied latitude and longitude data is in the same order as the distance data.')
  }


  if (length(RangeSamp)>2){
    stop('Error: the range sampling provided contains more than 2 values.
     \n Please provide either 1 value if you wish to sample both latitude and longitude equally
     \n (i.e. if the region you are looking at is approximately square) or
     \n please provide 2 values if you wish them to be sampled to different levels
     \n (i.e. if the region you are looking at is not square')
  }

  if (length(LongRange)!=2){
    stop('Error: you have provided more than 2 values in LongRange, these should be the 2 values denoting the maximum and minimum longitude range to be examined')

  }

  if (length(LatRange)!=2){
    stop('Error: you have provided more than 2 values in LatRange, these should be the 2 values denoting the maximum and minimum latitude range to be examined')
  }


  if (!(is.character(Method) | length(Method)>1)){

    stop('Error: you have not selected your preferred correlation Method.
     \n Please provide a single method, either Spearman or Pearson')

  }

}




