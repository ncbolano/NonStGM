#' Transformation of frequencies
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @return
#' @noRd
J.to.J.k.nu = function(J, k, nu) {
  n = nrow(J)
  J1 = J[(k + (-nu:nu)), ]
  J2 = c(t(J1))
  return(J2)
}

#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
J.to.J.k.nu.a = function(J, k, nu, a) {
  n = nrow(J)
  p = ncol(J)
  Jtemp = c(t(J[(k + (-nu:nu)), ]))
  Jtemp2 = Jtemp[-(p * nu + a)]
  return(Jtemp2)
}
