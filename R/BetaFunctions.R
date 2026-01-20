#' Transformation of frequencies
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param W weight Kernel chosen from (Uniform , Triangular, Quadratic)
#' @param M scalar smoothing value
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
beta = function(J, k, nu,W, M, a, delta) {
  HatR = Hat.R.reduced(J, k, nu, W,M, a)
  if(delta>0)
  {
    dim1 = ncol(HatR)
    HatR = HatR+ diag(rep(delta, dim1))

  }
  hatr = Hat.r(J, k, nu, W,M, a)
  betaC = solve(a = HatR,b = hatr)
  return(beta_beta_r(betaC,nu,ncol(J)))
}

beta_beta_r = function(beta,nu,p)
{
  tmp1 = as.vector(unlist(sapply((-nu):nu,function(i)
  {
    if(i==0)
      return(
        (nu)*p+(1:(p-1))
      )
    if(i<0)
      return(
        (abs(i)+nu)*p+1:p-1
      )
    if(i>0)
      return(
        (nu-abs(i))*p+1:p
      )
  }
  )))
  tmp2 = seq_along(beta)
  beta_erg = beta
  beta_erg[tmp2!=tmp1] = ((beta+Conj(beta[tmp1]))[tmp2!=tmp1])/2
  return(beta_erg)
}

beta_cv2 = function(J, k, nu,W, M_grid, a, delta)
{
  n_points = floor(median(M_grid)/(2*nu+1))
  tmp = sapply(M_grid,function(M)
  {
    J_tmp = J
    #Computational less demanding; one beta for all left out frequencies
    #more demanding, each left ouf beta gets its own frequency
    Eval= sapply(k+(-n_points):n_points*(2*nu+1),function(k_i)
    {
      J_tmp[k_i,a] = 0
      HatR = Hat.R.reduced(J, k_i, nu,W, M, a)
      hatr = Hat.r(J_tmp, k_i, nu,W, M, a)
      beta_tmp = solve(a = HatR,b = hatr)
      J[k_i, a]-J.to.J.k.nu.a(J, k_i, nu, a)%*%Conj(beta_tmp)
    })
    return(mean(abs(Eval)^2))
  })
  names(tmp)=M_grid
  return(tmp)
}

