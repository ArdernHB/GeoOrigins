
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

chr2nu <- function(X){as.numeric(as.chracter(X))}








