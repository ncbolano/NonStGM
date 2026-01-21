#' #' Transformation of frequencies
#' #'
#' #' @param J P-dimensional discrete fourier transform
#' #' @param k scalar value for \omega(k) = (2pi *k*) / n
#' #' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' #' @return
#' #' @noRd
#' J.to.J.k.nu = function(J, k, nu) {
#'   n = nrow(J)
#'   J1 = J[(k + (-nu:nu)), ]
#'   J2 = c(t(J1))
#'   return(J2)
#' }

#' #' Transformation of frequencies (excluding a single index)
#' #'
#' #' @param J P-dimensional discrete fourier transform
#' #' @param k scalar value for \omega(k) = (2pi *k*) / n
#' #' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' #' @param a scalar index value which determines removal of specific index from local DFT matrix
#' #' @return
#' #' @noRd
#' J.to.J.k.nu.a = function(J, k, nu, a) {
#'   n = nrow(J)
#'   p = ncol(J)
#'   Jtemp = c(t(J[(k + (-nu:nu)), ]))
#'   Jtemp2 = Jtemp[-(p * nu + a)]
#'   return(Jtemp2)
#' }


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

  indices_to_get = k + (-nu:nu)

  #    This maps any index outside of 1:n back into the valid range.
  wrapped_indices = ((indices_to_get - 1) %% n) + 1

  Jtemp = c(t(J[wrapped_indices, ]))

  Jtemp2 = Jtemp[-(p * nu + a)]
  return(Jtemp2)
}


#' #' Transformation of frequencies
#' #'
#' #' @param J P-dimensional discrete fourier transform
#' #' @param k scalar value for \omega(k) = (2pi *k*) / n
#' #' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' #' @return
#' #' @noRd
J.to.J.k.nu = function(J, k, nu) {
  n = nrow(J)

  indices_to_get = k + (-nu:nu)
  wrapped_indices = ((indices_to_get - 1) %% n) + 1

  J1 = J[wrapped_indices, ]
  J2 = c(t(J1))
  return(J2)
}
