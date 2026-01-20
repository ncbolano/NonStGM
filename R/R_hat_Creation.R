#' R_hat matrix creation associated with J.to.J.k.nu
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param W weight Kernel chosen from (Uniform , Triangular, Quadratic)
#' @param M scalar smoothing value
#' @return
#' @noRd
Hat.R = function(J, k, nu, W, M) {
  n = nrow(J)
  p = ncol(J)

  R_nu_M = matrix(0, p * (2 * nu + 1), p * (2 * nu + 1))

  for (s in (-M:M))
  {
    J_k_S = J.to.J.k.nu(J, k + s, nu)
    R_nu_M = R_nu_M + W(s/M)*J_k_S %*% t(Conj(J_k_S))
  }
  norm=sum(sapply((-M:M),function(s) W(s/M)))
  return(R_nu_M / norm)
}

#' R_hat matrix creating associated with J.to.J.k.nu.a
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param W weight Kernel chosen from (Uniform , Triangular, Quadratic)
#' @param M1 scalar smoothing value
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
Hat.R.reduced = function(J, k, nu, W,M1, a) {
  n = nrow(J)
  p = ncol(J)
  R_nu_M = matrix(0, (p * (2 * nu + 1) - 1), (p * (2 * nu + 1) - 1))

  for (s in (-M1:M1))
  {
    J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
    R_nu_M = R_nu_M + W(s/M1)*J_k_S %*% t(Conj(J_k_S))
  }
  norm=sum(sapply((-M1:M1),function(s) W(s/M1)))
  return(R_nu_M / norm)
}

#' R_hat vector creation excluding single missing frequency as computed in J.to.J.k.nu.a

#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param W weight Kernel chosen from (Uniform , Triangular, Quadratic)
#' @param M scalar smoothing value
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
Hat.r = function(J, k, nu,W, M, a) {
  n = nrow(J)
  p = ncol(J)

  r_nu_M = rep(0, (p * (2 * nu + 1) - 1))

  for (s in (-M:M))
  {
    J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
    r_nu_M = r_nu_M + W(s/M)* J_k_S * Conj(J[k + s, a])
  }
  norm=sum(sapply((-M:M),function(s) W(s/M)))
  return(r_nu_M / norm)
}

