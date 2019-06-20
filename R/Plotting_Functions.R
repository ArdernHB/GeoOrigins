
#' Internal function: Transparent named colour
#'
#' This function takes a named colour and returns the transparent equivalent
#' @param Colour A colour name from colours() function which is desired in transparent form.
#' @param alpha The level of transparency from 1 (completely transparent) to 100 (completely opaque) that the returned colour should be.
#' @return The transparent equivalent of a named colour
#' @keywords internal
#' @keywords colour
#' @keywords transparency
#' @author Ardern Hulme-Beaman
#' @examples
#' transpar(Colour = 'red', alpha = 50)

transpar<-function(Colour, alpha=100){
  newColour<-col2rgb(Colour)
  apply(newColour, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


#' Internal function: Constructing smooth contour polygon
#'
#' This function takes the Lat Longs that form a range around which a smooth contour is required. The contour is calculated base on kernal density estimation. This function uses the \code{kde} function of the \code{ks} Kernel Smoothing package and \code{contourLines} function of \code{grDevices} package.
#' @param LatLongs A dataframe with a 'Lat' column and a 'Long' column containing all latitude-longitude co-ordinates that are to form the basis upon which the contour is formed.
#' @return Contour polygon
#' @author Anna Rudzinski
#' @author Ardern Hulme-Beaman
#' @keywords internal
#' @keywords contour

Construct_contour <- function(LatLongs) {
  x <-  cbind(chr2nu(LatLongs$Longs), chr2nu(LatLongs$Lats))
  KSres <- ks::kde(x=x, H=Hpi(x=x), compute.cont=TRUE)
  contour.95 <- with(KSres, grDevices::contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate,levels=cont["5%"])[[1]])
  return(contour.95)
}


