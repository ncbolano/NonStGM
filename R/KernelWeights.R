#' Usage of a Uniform kernel for weighting coefficients
#'
#' @param x Input
#'
#' @return Returns 1 for values ( in absolute terms) less or equal 1, else 0.
#' @noRd
Kernel_Uniform = function(x) {
  return(as.numeric(abs(x)<=1))
  }



#' Usage of Triangular kernel for weighting coefficients
#'
#' @param x Input index
#'
#' @return Returns 1-|x| for values ( in absolute terms) less or equal 1, else 0.
#' @noRd
Kernel_Triangular = function(x) {
  return(as.numeric(abs(x)<=1)*(1-abs(x)))
}

#' Usage of Quadratic kernel for weighting coefficients
#'
#' @param x Input
#'
#' @return Returns 1-(|x|^2) for values ( in absolute terms) less or equal 1, else 0.
#' @noRd
Kernel_Quadtratic=function(x)
{
  return(as.numeric(abs(x)<=1)*(1-abs(x)^2))
}
