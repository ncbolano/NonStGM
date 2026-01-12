#' Usage of a Uniform kernel for weighting coefficients
#'
#' @param x Input index
#'
#' @return Returns 1 for values ( in absolute terms) less or equal 1, else 0.
#' @export
#' @noRd
Kernel_Uniform = function(x) {
  return(as.numeric(abs(x)<=1))
  }



#' Usage of Triangular kernel for weighting coefficients
#'
#' @param x Input index
#'
#' @return Returns 1-|x| for values ( in absolute terms) less or equal 1, else 0.
#' @export
#' @noRd
Kernel_Triangular = function(x) {
  return(as.numeric(abs(x)<=1)*(1-abs(x)))
}

