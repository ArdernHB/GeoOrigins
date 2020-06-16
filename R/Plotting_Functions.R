
#' Internal function: Transparent named colour
#'
#' This function takes a named colour and returns the transparent equivalent
#' @param Colour A colour name from colours() function which is desired in transparent form.
#' @param alpha The level of transparency from 1 (completely transparent) to 100 (completely opaque) that the returned colour should be.
#' @return The transparent equivalent of a named colour
#' @keywords internal
#' @author Ardern Hulme-Beaman


transpar<-function(Colour, alpha=100){
  newColour<-grDevices::col2rgb(Colour)
  apply(newColour, 2, function(curcoldata){grDevices::rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


#' Internal function: Constructing smooth contour polygon
#'
#' This function takes the Lat Longs that form a range around which a smooth contour is required. The contour is calculated base on kernal density estimation. This function uses the \code{\link[ks]{kde}} function of the \code{ks} Kernel Smoothing package and \code{\link[grDevices]{contourLines}} function of \code{grDevices} package.
#' @param LatLongs A dataframe with a 'Lats' column and a 'Longs' column containing all latitude-longitude co-ordinates that are to form the basis upon which the contour is formed.
#' @return Contour polygon
#' @author Anna Rudzinski
#' @author Ardern Hulme-Beaman
#' @keywords internal
#' @import ks

Construct_contour <- function(LatLongs) {

  x <-  cbind(chr2nu(LatLongs$Longs), chr2nu(LatLongs$Lats))
  KSres <- ks::kde(x=x, H=ks::Hpi(x=x), compute.cont=TRUE)
  contour.95 <- with(KSres, grDevices::contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate,levels=cont["5%"])[[1]])
  return(contour.95)
}


#' Internal function: Comparison of permutation array with true values
#'
#' This function takes the results of the permutation correlation test and returns a matrix of logical values of whether
#' @param PermArray an array of the permutation correlation results where rows and columns correspond with Lats and Longs respectively and slices are permutation iterations.
#' @param TrueMatrix the matrix of correlation values calculated from the true values.
#' @param Sig a numeric value for the signficance level the t.test is to perform to.
#' @param Alt a character string to be passed to the alternative argument of the t.test.
#' @return Matrix of logical values for grid cell locations where the significance level is met.
#' @author Ardern Hulme-Beaman
#' @keywords internal

PermArraySig <- function(PermArray, TrueMatrix, Sig, Alt='two.sided') {
  #Matrix=TrueVarValsMat; Array=PermVarHeat; Sig=SigLevel

  AdjMatrix <- matrix(NA, nrow = nrow(TrueMatrix), ncol = ncol(TrueMatrix))

  for (i in 1:dim(TrueMatrix)[1]){
    for(j in 1:dim(TrueMatrix)[2]){
      #i <- j <- 1
      TRes <- stats::t.test(x = PermArray[i,j,], mu = TrueMatrix[i,j], alternative = Alt)

      AdjMatrix[i,j] <- TRes$p.value<Sig
    }
  }

  return(AdjMatrix)
}

