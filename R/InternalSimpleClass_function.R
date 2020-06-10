
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
#' @keywords internal
#'
#' @author Ardern Hulme-Beaman




UserInputAssessment <- function(LatLongs, RefDistMat='skip', RefData='skip', DistVec='skip', Method){


  if (!(length(RefDistMat)==1)){
    if (dim(LatLongs)[1]!=dim(RefDistMat)[1]){
      stop('Error: the number of latitude and longitude coordinates you have provided do not match the number of specimens in the dataset.
     \n Please check these match and most importantly are in the same order and rerun the function')
    }
    if (dim(RefDistMat)[1]!=dim(RefDistMat)[2]){
      stop('Error: The distance matrix is not a square matrix. The distance matrix should be a square matrix of pairwise distances.
     \n Please check these match and most importantly are in the same order and rerun the function')
    }
  }

  if (!(length(DistVec)==1)){
    if (dim(LatLongs)[1]!=length(DistVec)){
      stop('Error: the number of latitude and longitude coordinates you have provided do not match the number of distances provided.
     \n Please check these match and most importantly are in the same order and rerun the function')
    }
  }



  if (!(length(RefData)==1)){
    if (length(dim(RefData))==2){
      specimenNo <- 1
    } else {
      specimenNo <- 3
    }
    if (dim(LatLongs)[1]!=dim(RefData)[specimenNo]){
      stop('Error: the number of latitude and longitude coordinates you have provided do not match the number of specimens in the dataset.
     \n Please check these match and most importantly are in the same order and rerun the function')
    }
  }



  if (dim(LatLongs)[2]!=2){
    stop('Error: the matrix of latitude and longitude coordinates you have provided does not contain 2 columns.
     \n Please check that you have supplied the correct data and in the correct format.
     \n Importantly, please check that the supplied latitude and longitude data is in the same order as the distance data.')
  }



  if (!(is.character(Method) | length(Method)>1)){

    stop('Error: you have not selected your preferred correlation Method.
     \n Please provide a single method, either Spearman or Pearson')

  }

}




