################# extractVariables.R

#' Extracts significant variables from original series
#'
#' @param x multivariate time series with dimension nxp
#' @return list of three variables (JJ,p,nu)
#' @noRd
extractVariables = function(x) {
  p = ncol(x)
  nu = 2
  JJ0 = mvfft(x) / sqrt(nrow(x))
  n_rows = nrow(x)
  backward = c(((n_rows / 2) + 1):n_rows)
  JJ = rbind(
    JJ0[backward, , drop = FALSE],
    JJ0[1:(n_rows / 2), , drop = FALSE]
  )
  return(list(JJ,p,nu))
}

########## DFTransform.R
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






################## Combined_MK_Estimation.R


#' Assigning r values to a specific index location
#'
#' @param r P-dimensional discrete fourier transform
#' @param c scalar value for \omega(k) = (2pi *k*) / n
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param p Dimension multivariate time series (# nodes)
#' @return
#' @noRd
r_to_loc = function(r, c, a, nu, p) {
  if (r < 0) {
    return(p * (r + nu) + c)
  }
  if (r == 0 & c > a) {
    return(p * nu + (c - 1))
  }
  if (r == 0 & c < a) {
    return(p * nu + c)
  }
  if (r > 0) {
    return(p * nu + (p - 1) + (r - 1) * p + c)
  }
}

#' Transformation of frequencies (excluding a single index)
#'
#' @param S_hat
#' @param I_k Identity matrix of square dimension k
#' @return
#' @noRd
#'
whittle_loss = function(S_hat, I_k) {
  S_hat = Re(S_hat)

  # Regularization for invertibility if needed
  S_hat = S_hat + diag(1e-6, nrow(S_hat))

  return( Re(sum(diag(solve(S_hat) %*% I_k)) + log(det(S_hat))))
}

#' Transformation of frequencies (excluding a single index)
#'
#' @param JJ P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param M scalar value for smoothing distance (full smoothing dist. 2M)
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
smoothed_spectral_density_LOO = function(JJ, k, M, Kernel_func = Kernel_Triangular, leave_out = k) {
  n = nrow(JJ)
  p = ncol(JJ)
  S_k = matrix(0+0i, p, p)
  weight_sum = 0

  for (j_offset in -M:M) {
    current_j = k + j_offset

    idx = ((current_j - 1) %% n) + 1

    if (idx == leave_out) next

    weight = Kernel_func(j_offset / M)

    if (weight > 0) {
      I_j = outer(JJ[idx, ], Conj(JJ[idx, ]))
      S_k = S_k + weight * I_j
      weight_sum = weight_sum + weight
    }
  }

  if (weight_sum > 0) S_k = S_k / weight_sum
  return(S_k)
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
smoothed_spectral_density = function(JJ, k, M, Kernel_func = Kernel_Triangular) {
  n = nrow(JJ)
  p = ncol(JJ)
  S_k = matrix(0 + 0i, p, p)
  weight_sum = 0

  for (j_offset in -M:M) {
    current_j = k + j_offset

    idx = ((current_j - 1) %% n) + 1

    weight = Kernel_func(j_offset / M)

    if (weight > 0) {
      I_j = outer(JJ[idx, ], Conj(JJ[idx, ]))
      S_k = S_k + weight * I_j
      weight_sum = weight_sum + weight
    }
  }

  if (weight_sum > 0) {
    S_k = S_k / weight_sum
  }

  return(S_k)
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
return_max = function(x, M) {
  max1 = which.max(x)

  # Create a copy and set values within 2*M of max1 to -Inf
  x2 = x
  min_dist =  M

  # Set forbidden region to -Inf
  forbidden_indices = max(1, max1 - min_dist):min(length(x), max1 + min_dist)
  x2[forbidden_indices] = -Inf

  # Find second maximum (at least 2*M away)
  max2 = which.max(x2)

  # Validity check of existing other point to prevent breakdown for M > n/4
  if (is.infinite(x2[max2])) {
    warning("No valid second maximum found at distance >= 2*M")
    return(c(max1, NA))
  }

  return(c(max1, max2))
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
extractK = function(JJ,coefnum) {
  cat("\n--- extractK ---\n")
  cat("dim(JJ):", paste(dim(JJ), collapse=" x "), "\n")
  cat("coefnum:", coefnum, "\n")
  M_initial = 2 * coefnum
  n_freq = floor(nrow(JJ) / 2)
  floor = M_initial + 1
  ceiling = n_freq - M_initial
  cat("M_initial =", M_initial, "\n")
  cat("n_freq =", n_freq, "\n")
  cat("floor =", floor, "\n")
  cat("ceiling =", ceiling, "\n")

  trace_smooth = numeric(n_freq)
  diag_smooth = matrix(0, nrow = n_freq, ncol = ncol(JJ)) # Matrix for diagonal elements
  largest_eig_smooth = numeric(n_freq)
  condition_number_smooth = numeric(n_freq)

  for (k in floor:ceiling) {

    # Smoothed Periodogram (weighted by kernel)
    S_k = Re(smoothed_spectral_density(JJ, k, M_initial, Kernel_Triangular))
    trace_smooth[k] = sum(diag(S_k))
    diag_smooth[k, ] = diag(S_k)

    eigenvals_smooth = eigen(S_k, only.values = TRUE)$values
    eigenvals_smooth = Re(eigenvals_smooth)
    largest_eig_smooth[k] = max(eigenvals_smooth)
    condition_number_smooth[k] = max(eigenvals_smooth) / (min(eigenvals_smooth) + 1e-10)
  }

  composite_smooth = trace_smooth / ((condition_number_smooth)^(.5)) # From frequency 1:(n/2)
  chosen_frequencies = return_max(composite_smooth , M_initial)
  k = chosen_frequencies
  return(k)
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
local_M_selection = function(JJ, k, M_grid) {

  best_local_M = numeric(length = length(k))
  CV_scores = numeric(length = length(M_grid))
  for (i in seq_along(k)) {
    k_candidate = k[i]
    I_candidate_k = outer(JJ[k_candidate, ], Conj(JJ[k_candidate, ]))

    for (j in seq_along(M_grid)) {
      M = M_grid[j]

      S_hat_k = smoothed_spectral_density_LOO(JJ, k = k_candidate, M = M, leave_out = k_candidate)
      CV_scores[j] = Re(whittle_loss(S_hat_k, I_candidate_k))

      #cat("M =", M, "K =", k_candidate , " CV =", CV_scores[j], "\n")
    }

    best_M_for_this_k = M_grid[which.min(CV_scores)]

    best_local_M[i] = best_M_for_this_k
  }
  print(best_local_M)
  return(best_local_M)
}

#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
extractM = function(JJ,k,coefnum) {
  n = nrow(JJ)
  M_grid = unique(round(seq(max(coefnum , n^(1/5)), max(coefnum , n^(1/5)) + 1, length.out = 10)))
  M_list = local_M_selection(JJ, k, M_grid)

  return(M_list)
}

############################# KernelWeights.R

#' Usage of a Uniform kernel for weighting coefficients
#'
#' @param x Input
#'
#' @return Returns 1 for values ( in absolute terms) less or equal 1, else 0.
#' @noRd
Kernel_Uniform = function(x) {
  return(as.numeric(abs(x)<=1))
}






############################# R_hat_Creation.R

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
Kernel_Quadratic=function(x)
{
  return(as.numeric(abs(x)<=1)*(1-abs(x)^2))
}






####################### BetaFunctions.R

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
beta = function(J, k, nu, W, M, a, delta) {
  HatR = Hat.R.reduced(J, k, nu, W, M, a)
  if(delta>0)
  {
    dim1 = ncol(HatR)
    HatR = HatR+ diag(rep(delta, dim1))

  }
  # Added check
  if (rcond(HatR) < .Machine$double.eps) {
    p = ncol(J)
    len_out = p * (2 * nu + 1) - 1
    return(rep(NA, len_out))
  }

  hatr = Hat.r(J, k, nu, W,M, a)
  betaC = solve(a = HatR,b = hatr)
  return(beta_beta_r(betaC,nu,ncol(J)))
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
beta_beta_r = function(beta, nu, p)
{
  tmp1 = as.vector(unlist(sapply((-nu):nu, function(i)
  {
    if (i == 0)
      return((nu) * p + seq_len(p - 1))
    if (i < 0)
      return((abs(i) + nu) * p + 1:p - 1)
    if (i > 0)
      return((nu - abs(i)) * p + 1:p)
  })))
  tmp2 = seq_along(beta)
  beta_erg = beta
  beta_erg[tmp2 != tmp1] = ((beta + Conj(beta[tmp1]))[tmp2 != tmp1]) / 2
  return(beta_erg)
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
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





######## VarianceFunctions.R
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
variance.estimator.Re.v2 = function(J, k, nu, M, a, L, delta) {
  betaC = beta(J, k, nu, M, a, delta)
  dim1 = length(betaC)
  HatRred = Hat.R.reduced(J, k, nu, M, a)
  ridge = diag(rep(delta, dim1))
  invHatRreduced = solve(HatRred + ridge)
  # invHatRreduced = solve(Hat.R.reduced(J,k,nu,M1=M,a)+diag(rep(delta,dim1)))
  #          dim1 = ncol(invHatRreduced)
  Y = matrix(rep(0, (2 * M + 1) * dim1), nrow = dim1)
  p = ncol(J)

  for (s in (-M:M))
  {
    J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
    temp1 = (Conj(J[k + s, a]) - sum(betaC * Conj(J_k_S)))
    temp2 = invHatRreduced %*% J_k_S
    Y[, s + M + 1] = Re(temp1 * temp2 / (2 * M + 1))
  }

  sigma1 = matrix(rep(0, dim1**2), ncol = dim1)

  for (s1 in (1:(2 * M + 1))) {
    for (s2 in (1:(2 * M + 1))) {
      if (abs(s1 - s2) < L + 1) sigma1 = sigma1 + Y[, s1] %*% t(Y[, s2])
    }
  }

  sigma2 = matrix(rep(0, dim1**2), ncol = dim1)
  k1 = (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 = sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }

  sigmaRe = diag(sigma1 + sigma2)

  return(sigmaRe)
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
variance.estimator.Im.v2 = function(J, k, nu, M, a, L, delta) {
  betaC = beta(J, k, nu, M, a, delta)
  dim1 = length(betaC)
  HatRred = Hat.R.reduced(J, k, nu, M, a)
  ridge = diag(rep(delta, dim1))
  invHatRreduced = solve(HatRred + ridge)
  # invHatRreduced = solve(Hat.R.reduced(J,k,nu,M1=M,a)+diag(rep(delta,dim1)))
  #          dim1 = ncol(invHatRreduced)
  Y = matrix(rep(0, (2 * M + 1) * dim1), nrow = dim1)
  p = ncol(J)

  for (s in (-M:M))
  {
    J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
    temp1 = (Conj(J[k + s, a]) - sum(betaC * Conj(J_k_S)))
    temp2 = invHatRreduced %*% J_k_S
    Y[, s + M + 1] = Im(temp1 * temp2 / (2 * M + 1))
  }

  sigma1 = matrix(rep(0, dim1**2), ncol = dim1)

  for (s1 in (1:(2 * M + 1))) {
    for (s2 in (1:(2 * M + 1))) {
      if (abs(s1 - s2) < L+1) sigma1 = sigma1 + Y[, s1] %*% t(Y[, s2])
    }
  }

  sigma2 = matrix(rep(0, dim1**2), ncol = dim1)
  n = nrow(J)
  k1 = (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 = sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }

  sigmaIm = diag(sigma1 + sigma2)

  return(sigmaIm)
}

#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
variance.estimator.v2 = function(J, k, nu,W, M, a, L, delta) {
  betaC = beta(J, k, nu, W,M, a, delta)
  dim1 = length(betaC)
  HatRred = Hat.R.reduced(J, k, nu,W, M, a)
  ridge = diag(rep(delta, dim1))
  invHatRreduced = solve(HatRred + ridge)
  norm=sum(sapply((-M:M),function(s) W(s/M)))
  # invHatRreduced = solve(Hat.R.reduced(J,k,nu,M1=M,a)+diag(rep(delta,dim1)))
  #          dim1 = ncol(invHatRreduced)
  Y = matrix(rep(0, 2*(2 * M + 1) * dim1), nrow = 2*dim1)
  p = ncol(J)

  for (s in (-M:M))
  {
    J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
    temp1 = (Conj(J[k + s, a]) - sum(betaC * Conj(J_k_S)))
    temp2 = invHatRreduced %*% J_k_S
    Y[, s + M + 1] = W(s/M)*c(Re(temp1 * temp2 / norm),Im(temp1 * temp2 / norm))
  }

  sigma1 = matrix(rep(0, (2*dim1)**2), ncol = 2*dim1)

  for (s1 in (1:(2 * M + 1))) {
    for (s2 in (1:(2 * M + 1))) {
      if (abs(s1 - s2) < L + 1) sigma1 = sigma1 + Y[, s1] %*% t(Y[, s2])
    }
  }

  sigma2 = matrix(rep(0, (2*dim1)**2), ncol = 2*dim1)
  n = nrow(J)
  k1 = (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 = sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }

  return(sigma1+sigma2)
}




################## parallel_sim_dependencies.R


# ============================================================
# p=3 VAR SIMULATION FUNCTION
# ============================================================

library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)

p = 3

sim.tvVAR_p3 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5, 0.2, 0,
    0,   0.8, 0,
    0,   0.3, 0.6
  ), ncol = 3, byrow = TRUE)

  n = m + burnin
  p = 3
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)

  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }

  x2 = x1[-c(1:burnin), ]
  return(x2)
}



check_tvVAR_stability_p3 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5, 0.2, 0,
    0,   0.8, 0,
    0,   0.3, 0.6
  ), ncol = 3, byrow = TRUE)

  n = burnin + m
  st = 0.3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  max_mod = sapply(st, function(s) { A[1, 1] = s; max(Mod(eigen(A)$values)) })
  cat("A[1,1] range:", round(range(st), 3), "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}

# ============================================================
# GRAPH RECOVERY FUNCTION FOR p=3
# ============================================================
# True edges:
#   Time-varying: (1,1)
#   Off-diagonal: (1,2), (2,3)
#   Total: 3 true edges, 2 off-diagonal

graph_recovery_p3 = function(Test_tibble, alpha = 0.05) {

  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r != 0))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  graph_recovery_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      # True edges
      detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
      detect_12 = any(a == 1 & c == 2 & p_min < alpha),
      detect_23 = any(a == 2 & c == 3 & p_min < alpha),
      # False edges (null)
      detect_13 = any(a == 1 & c == 3 & p_min < alpha),
      detect_22 = any(a == 2 & c == 2 & p_min < alpha),
      detect_33 = any(a == 3 & c == 3 & p_min < alpha),
      # Edge counts
      n_offdiag_edges = sum(a != c & p_min < alpha),
      n_diag_edges = sum(a == c & p_min < alpha),
      .groups = "drop"
    )

  accuracy_summary = graph_recovery_results |>
    mutate(
      all_offdiag_detected = detect_12 & detect_23,
      all_true_detected = all_offdiag_detected & detect_11_tv,
      perfect = all_offdiag_detected & (n_offdiag_edges == 2),
      sensitivity = (detect_11_tv + detect_12 + detect_23) / 3,
      false_positives = detect_13 + detect_22 + detect_33
    )

  overall_accuracy = tibble(
    Method = "Global BY",
    Perfect_Recovery = mean(accuracy_summary$perfect),
    All_True_Detected = mean(accuracy_summary$all_true_detected),
    # True edges
    Detect_11_TV = mean(accuracy_summary$detect_11_tv),
    Detect_12 = mean(accuracy_summary$detect_12),
    Detect_23 = mean(accuracy_summary$detect_23),
    # False edges
    False_13 = mean(accuracy_summary$detect_13),
    False_22 = mean(accuracy_summary$detect_22),
    False_33 = mean(accuracy_summary$detect_33),
    # Summary
    Mean_Sensitivity = mean(accuracy_summary$sensitivity),
    Mean_N_Offdiag_Edges = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Positives = mean(accuracy_summary$false_positives)
  )

  return(overall_accuracy)
}


#p = 5

sim.tvVAR_p5 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,
    0.3,  0.6,  0,    0,    0,
    0,    0.3,  0.5,  0,    0,
    0,    0,    0.3, 0.55, 0,
    0.3, 0,    0,    0.3,  0.5
  ), ncol = 5, byrow = TRUE)
  n = m + burnin
  p = 5
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)
  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}

check_tvVAR_stability = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,
    0.3,  0.6,  0,    0,    0,
    0,    0.3,  0.5,  0,    0,
    0,    0,    0.3, 0.55, 0,
    0.3, 0,    0,    0.3,  0.5
  ), ncol = 5, byrow = TRUE)
  n = burnin + m
  st = 0.3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  max_mod = sapply(st, function(s) { A[1,1] = s; max(Mod(eigen(A)$values)) })
  cat("A[1,1] range:", round(range(st), 3), "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}

# graph_recovery_p5 = function(Test_tibble, alpha = 0.05) {
#   Test_tibble = Test_tibble |>
#     filter(a <= c) |>
#     filter((a != c & r == 0) | (a == c & r != 0))
#
#   Graph_tibble = Test_tibble |>
#     group_by(a, c, i) |>
#     summarise(
#       p_min = min(p.adjust(c(Re, Im), method = "BY")),
#       .groups = "drop"
#     )
#
#   graph_recovery_results = Graph_tibble |>
#     group_by(i) |>
#     summarise(
#       detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
#       detect_12 = any(a == 1 & c == 2 & p_min < alpha),
#       detect_14 = any(a == 1 & c == 4 & p_min < alpha),
#       detect_15 = any(a == 1 & c == 5 & p_min < alpha),
#       detect_23 = any(a == 2 & c == 3 & p_min < alpha),
#       detect_34 = any(a == 3 & c == 4 & p_min < alpha),
#       detect_45 = any(a == 4 & c == 5 & p_min < alpha),
#       n_offdiag_edges = sum(a != c & p_min < alpha),
#       .groups = "drop"
#     )
#
#   accuracy_summary = graph_recovery_results |>
#     mutate(
#       all_offdiag_detected = detect_12 & detect_14 & detect_15 &
#         detect_23 & detect_34 & detect_45,
#       all_true_detected = all_offdiag_detected & detect_11_tv,
#       perfect = all_true_detected & (n_offdiag_edges == 6),
#       sensitivity = (detect_11_tv + detect_12 + detect_14 +
#                        detect_15 + detect_23 + detect_34 + detect_45) / 7,
#       false_positives = pmax(0, n_offdiag_edges - 6)
#     )
#
#   overall_accuracy = tibble(
#     Method = "Global BY",
#     Perfect_Recovery = mean(accuracy_summary$perfect),
#     All_True_Detected = mean(accuracy_summary$all_true_detected),
#     Detect_11_TV = mean(accuracy_summary$detect_11_tv),
#     Detect_12 = mean(accuracy_summary$detect_12),
#     Detect_14 = mean(accuracy_summary$detect_14),
#     Detect_15 = mean(accuracy_summary$detect_15),
#     Detect_23 = mean(accuracy_summary$detect_23),
#     Detect_34 = mean(accuracy_summary$detect_34),
#     Detect_45 = mean(accuracy_summary$detect_45),
#     Mean_Sensitivity = mean(accuracy_summary$sensitivity),
#     Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
#     Mean_False_Positives = mean(accuracy_summary$false_positives)
#   )
#   return(overall_accuracy)
# }
# ============================================================
# SIMULATION RUNNER FOR p=5
# ============================================================
run_NonStGM_simulation_p5 = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                     TV_size = 0.6, nu = 2, L = 1,
                                     Kernel = 'Kernel_Triangular', seed = 0) {
  All_Test_tibble = NULL
  set.seed(seed)
  for (iter in 1:R) {
    x = sim.tvVAR_p5(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    All_Test_tibble = rbind(All_Test_tibble, tib)
  }
  results = graph_recovery_p5(All_Test_tibble, alpha = alpha)
  return(results)
}

# results_p5 = run_NonStGM_simulation_p5(
#   R = 100,
#   alpha = 0.05,
#   burnin = 200,
#   m = 2^12,
#   TV_size = 0.6,
#   nu = 2,
#   L = 1,
#   Kernel = 'Kernel_Triangular',
#   seed = 1
# )
# print(results_p5)

# ============================================================
# p=7 VAR SIMULATION FUNCTION
# ============================================================
p = 7

sim.tvVAR_p7 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,    0,    0,
    0.3,  0.5,  0,    0,    0,    0,    0,
    0,    0.3,  0.5,  0,    0,    0,    0,
    0,    0,    0.3,  0.5,  0,    0,    0,
    0,    0,    0,    0.3,  0.5,  0,    0,
    0,    0,    0,    0,    0.3,  0.5,  0,
    0.3,  0,    0,    0,    0,    0.3,  0.5
  ), ncol = 7, byrow = TRUE)
  n = m + burnin
  p = 7
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = .3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1) # CHANGED
  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}

check_tvVAR_stability_p7 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,    0,    0,
    0.3,  0.5,  0,    0,    0,    0,    0,
    0,    0.3,  0.5,  0,    0,    0,    0,
    0,    0,    0.3,  0.5,  0,    0,    0,
    0,    0,    0,    0.3,  0.5,  0,    0,
    0,    0,    0,    0,    0.3,  0.5,  0,
    0.3,  0,    0,    0,    0,    0.3,  0.5
  ), ncol = 7, byrow = TRUE)
  n = burnin + m
  st = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1) #CHANGED
  max_mod = sapply(st, function(s) { A[1, 1] = s; max(Mod(eigen(A)$values)) })
  cat("A[1,1] range:", round(range(st), 3), "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}


graph_recovery_p7 = function(Test_tibble, alpha = 0.05) {
  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r == 1))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  graph_recovery_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
      detect_12 = any(a == 1 & c == 2 & p_min < alpha),
      detect_16 = any(a == 1 & c == 6 & p_min < alpha),
      detect_17 = any(a == 1 & c == 7 & p_min < alpha),
      detect_23 = any(a == 2 & c == 3 & p_min < alpha),
      detect_34 = any(a == 3 & c == 4 & p_min < alpha),
      detect_45 = any(a == 4 & c == 5 & p_min < alpha),
      detect_56 = any(a == 5 & c == 6 & p_min < alpha),
      detect_67 = any(a == 6 & c == 7 & p_min < alpha),
      n_offdiag_edges = sum(a != c & p_min < alpha),
      .groups = "drop"
    )

  accuracy_summary = graph_recovery_results |>
    mutate(
      all_offdiag_detected = detect_12 & detect_16 & detect_17 &
        detect_23 & detect_34 & detect_45 & detect_56 & detect_67,
      all_true_detected = all_offdiag_detected & detect_11_tv,
      perfect = all_true_detected & (n_offdiag_edges == 8),
      sensitivity = (detect_11_tv + detect_12 + detect_16 + detect_17 +
                       detect_23 + detect_34 + detect_45 + detect_56 + detect_67) / 9,
      false_positives = pmax(0, n_offdiag_edges - 8)
    )

  overall_accuracy = tibble(
    Method = "Global BY",
    Perfect_Recovery = mean(accuracy_summary$perfect),
    All_True_Detected = mean(accuracy_summary$all_true_detected),
    Detect_11_TV = mean(accuracy_summary$detect_11_tv),
    Detect_12 = mean(accuracy_summary$detect_12),
    Detect_16 = mean(accuracy_summary$detect_16),
    Detect_17 = mean(accuracy_summary$detect_17),
    Detect_23 = mean(accuracy_summary$detect_23),
    Detect_34 = mean(accuracy_summary$detect_34),
    Detect_45 = mean(accuracy_summary$detect_45),
    Detect_56 = mean(accuracy_summary$detect_56),
    Detect_67 = mean(accuracy_summary$detect_67),
    Mean_Sensitivity = mean(accuracy_summary$sensitivity),
    Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Positives = mean(accuracy_summary$false_positives)
  )
  return(overall_accuracy)
}
# ============================================================
# SIMULATION RUNNER FOR p=7
# ============================================================
run_NonStGM_simulation_p7 = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                     TV_size = 0.6, nu = 2, L = 1,
                                     Kernel = 'Kernel_Triangular', seed = 0) {
  All_Test_tibble = NULL
  set.seed(seed)
  for (iter in 1:R) {
    x = sim.tvVAR_p7(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    All_Test_tibble = rbind(All_Test_tibble, tib)
  }
  results = graph_recovery_p7(All_Test_tibble, alpha = alpha)
  return(list(results, All_Test_tibble))
}

# results_p7 = run_NonStGM_simulation_p7(
#   R = 100,
#   alpha = 0.05,
#   burnin = 20,
#   m = 2^12,
#   TV_size = 0.6,
#   nu = 2,
#   L = 1,
#   Kernel = 'Kernel_Triangular',
#   seed = 0
# )
# print(results_p7[[1]])
# hi = results_p7[[2]]




p = 7

sim.tvVAR2_p7 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.4,  0,    0,    0,    0,    0,    0.25,
    0.3,  0.4,  0,    0,    0,    0,    0,
    0,    0,    0.4,  0.3,  0,    0,    0,
    0,    0,    0,    0.4,  0.3,  0,    0,
    0,    0,    0,    0,    0.4,  0.3, 0,
    0,    0.25, 0,    0,    0,    0.4,  0.3,
    0,    0,    0,    0,    0,    0,    0.4
  ), ncol = 7, byrow = TRUE)
  n = m + burnin
  p = 7
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st3 = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  st4 = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  for (tt in 2:n) {
    A.t = A
    A.t[3, 3] = st3[tt]
    A.t[4, 4] = st4[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}
check_tvVAR2_stability_p7 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.4,  0,    0,    0,    0,    0,    0.25,
    0.3,  0.4,  0,    0,    0,    0,    0,
    0,    0,    0.4,  0.3,  0,    0,    0,
    0,    0,    0,    0.4,  0.3,  0,    0,
    0,    0,    0,    0,    0.4,  0.25, 0,
    0,    0.25, 0,    0,    0,    0.4,  0.3,
    0,    0,    0,    0,    0,    0,    0.4
  ), ncol = 7, byrow = TRUE)
  n = burnin + m
  st = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  max_mod = sapply(st, function(s) {
    A[3, 3] = s
    A[4, 4] = s
    max(Mod(eigen(A)$values))
  })
  cat("A[3,3]/A[4,4] range:", round(range(st), 3),
      "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}


graph_recovery2_p7 = function(Test_tibble, alpha = 0.05) {
  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r == 1))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  graph_recovery_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      detect_33_tv = any(a == 3 & c == 3 & p_min < alpha),
      detect_44_tv = any(a == 4 & c == 4 & p_min < alpha),
      detect_12    = any(a == 1 & c == 2 & p_min < alpha),
      detect_17    = any(a == 1 & c == 7 & p_min < alpha),
      detect_26    = any(a == 2 & c == 6 & p_min < alpha),
      detect_27    = any(a == 2 & c == 7 & p_min < alpha),
      detect_34    = any(a == 3 & c == 4 & p_min < alpha),
      detect_45    = any(a == 4 & c == 5 & p_min < alpha),
      detect_56    = any(a == 5 & c == 6 & p_min < alpha),
      detect_67    = any(a == 6 & c == 7 & p_min < alpha),
      n_offdiag_edges = sum(a != c & p_min < alpha),
      .groups = "drop"
    )

  accuracy_summary = graph_recovery_results |>
    mutate(
      all_offdiag_detected = detect_12 & detect_17 & detect_26 & detect_27 &
        detect_34 & detect_45 & detect_56 & detect_67,
      all_true_detected = all_offdiag_detected & detect_33_tv & detect_44_tv,
      perfect = all_true_detected & (n_offdiag_edges == 8),
      sensitivity = (detect_33_tv + detect_44_tv + detect_12 + detect_17 +
                       detect_26 + detect_27 + detect_34 + detect_45 +
                       detect_56 + detect_67) / 10,
      false_positives = pmax(0, n_offdiag_edges - 8)
    )

  overall_accuracy = tibble(
    Method               = "Global BY",
    Perfect_Recovery     = mean(accuracy_summary$perfect),
    All_True_Detected    = mean(accuracy_summary$all_true_detected),
    Detect_33_TV         = mean(accuracy_summary$detect_33_tv),
    Detect_44_TV         = mean(accuracy_summary$detect_44_tv),
    Detect_12            = mean(accuracy_summary$detect_12),
    Detect_17            = mean(accuracy_summary$detect_17),
    Detect_26            = mean(accuracy_summary$detect_26),
    Detect_27            = mean(accuracy_summary$detect_27),
    Detect_34            = mean(accuracy_summary$detect_34),
    Detect_45            = mean(accuracy_summary$detect_45),
    Detect_56            = mean(accuracy_summary$detect_56),
    Detect_67            = mean(accuracy_summary$detect_67),
    Mean_Sensitivity     = mean(accuracy_summary$sensitivity),
    Mean_N_Edges         = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Positives = mean(accuracy_summary$false_positives)
  )
  return(overall_accuracy)
}


p = 15
sim.tvVAR_p15 = function(burnin, m, TV_size) {
  A = matrix(0, ncol = 15, nrow = 15)
  diag(A) = 0.4
  for (i in 1:14) {
    A[i + 1, i] = 0.3
  }
  A[14,1] = 0.3
  n = m + burnin
  p = 15
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st1  = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  st12 = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  for (tt in 2:n) {
    A.t = A
    A.t[1, 1]   = st1[tt]
    A.t[12, 12] = st12[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}

check_tvVAR_stability_p15 = function(burnin, m, TV_size) {
  A = matrix(0, ncol = 15, nrow = 15)
  diag(A) = 0.4
  for (i in 1:14) {
    A[i + 1, i] = 0.3
  }
  n = burnin + m
  st = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  max_mod = sapply(st, function(s) {
    A[1, 1]   = s
    A[12, 12] = s
    max(Mod(eigen(A)$values))
  })
  cat("A[1,1]/A[12,12] range:", round(range(st), 3),
      "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}

graph_recovery_p15 = function(Test_tibble, alpha = 0.05) {
  # True graph: chain 1-2-3-...-15 (14 off-diagonal edges)
  # TV on diagonal: nodes 1 and 12
  # A'A only has nonzero off-diagonal at (j, j+1) since columns only overlap at one row
  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r == 1))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  graph_recovery_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      detect_1_1_tv  = any(a == 1  & c == 1  & p_min < alpha),
      detect_12_12_tv = any(a == 12 & c == 12 & p_min < alpha),
      detect_1_2   = any(a == 1  & c == 2  & p_min < alpha),
      detect_1_13   = any(a == 1  & c == 13  & p_min < alpha),
      detect_1_14   = any(a == 1  & c == 14  & p_min < alpha),
      detect_2_3   = any(a == 2  & c == 3  & p_min < alpha),
      detect_3_4   = any(a == 3  & c == 4  & p_min < alpha),
      detect_4_5   = any(a == 4  & c == 5  & p_min < alpha),
      detect_5_6   = any(a == 5  & c == 6  & p_min < alpha),
      detect_6_7   = any(a == 6  & c == 7  & p_min < alpha),
      detect_7_8   = any(a == 7  & c == 8  & p_min < alpha),
      detect_8_9   = any(a == 8  & c == 9  & p_min < alpha),
      detect_9_10  = any(a == 9  & c == 10 & p_min < alpha),
      detect_10_11 = any(a == 10 & c == 11 & p_min < alpha),
      detect_11_12 = any(a == 11 & c == 12 & p_min < alpha),
      detect_12_13 = any(a == 12 & c == 13 & p_min < alpha),
      detect_13_14 = any(a == 13 & c == 14 & p_min < alpha),
      detect_14_15 = any(a == 14 & c == 15 & p_min < alpha),
      n_offdiag_edges = sum(a != c & p_min < alpha),
      .groups = "drop"
    )

  accuracy_summary = graph_recovery_results |>
    mutate(
      all_offdiag_detected = detect_1_2 & detect_1_13 & detect_1_14 & detect_2_3 & detect_3_4 &
        detect_4_5 & detect_5_6 & detect_6_7 & detect_7_8 &
        detect_8_9 & detect_9_10 & detect_10_11 & detect_11_12 &
        detect_12_13 & detect_13_14 & detect_14_15,
      all_true_detected = all_offdiag_detected & detect_1_1_tv & detect_12_12_tv,
      perfect = all_true_detected & (n_offdiag_edges == 14),
      sensitivity = (detect_1_1_tv + detect_12_12_tv +
                       detect_1_2 + detect_1_14 + detect_1_13 + detect_2_3 + detect_3_4 + detect_4_5 +
                       detect_5_6 + detect_6_7 + detect_7_8 + detect_8_9 +
                       detect_9_10 + detect_10_11 + detect_11_12 +
                       detect_12_13 + detect_13_14 + detect_14_15) / 18,
      false_positives = pmax(0, n_offdiag_edges - 16)
    )

  overall_accuracy = tibble(
    Method               = "Global BY",
    Perfect_Recovery     = mean(accuracy_summary$perfect),
    All_True_Detected    = mean(accuracy_summary$all_true_detected),
    Detect_1_1_TV        = mean(accuracy_summary$detect_1_1_tv),
    Detect_12_12_TV      = mean(accuracy_summary$detect_12_12_tv),
    Detect_1_2           = mean(accuracy_summary$detect_1_2),
    Detect_2_3           = mean(accuracy_summary$detect_2_3),
    Detect_3_4           = mean(accuracy_summary$detect_3_4),
    Detect_4_5           = mean(accuracy_summary$detect_4_5),
    Detect_5_6           = mean(accuracy_summary$detect_5_6),
    Detect_6_7           = mean(accuracy_summary$detect_6_7),
    Detect_7_8           = mean(accuracy_summary$detect_7_8),
    Detect_8_9           = mean(accuracy_summary$detect_8_9),
    Detect_9_10          = mean(accuracy_summary$detect_9_10),
    Detect_10_11         = mean(accuracy_summary$detect_10_11),
    Detect_11_12         = mean(accuracy_summary$detect_11_12),
    Detect_12_13         = mean(accuracy_summary$detect_12_13),
    Detect_13_14         = mean(accuracy_summary$detect_13_14),
    Detect_14_15         = mean(accuracy_summary$detect_14_15),
    Mean_Sensitivity     = mean(accuracy_summary$sensitivity),
    Mean_N_Edges         = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Positives = mean(accuracy_summary$false_positives)
  )
  return(overall_accuracy)
}

NonStGM_adj = function(x, Kernel = 'Kernel_Triangular', nu = 2, L = 1, alpha = 0.05,
                       skeleton_only = FALSE) {

  Kernel_function = get(Kernel)

  variable_list = extractVariables(x)
  J = variable_list[[1]] ; p = as.numeric(variable_list[2])
  coefnum = (2 * p * nu) + p - 1

  # Adaptive frequency and smoothing distance choice for both methods (albeit first step) done here
  k = extractK(J, coefnum)
  n = nrow(x)
  M_grid = unique(round(seq(max(coefnum , n^(1/5)), max(coefnum + 100 , n^(1/2)) , length.out = 10)))
  M_list = local_M_selection(J, k, M_grid)

  # Empty array loading for Test_tibble filling
  n_k = length(k)
  betaCoefAll = array(0, c(n_k, coefnum, p))
  varbeta = array(0, c(n_k, coefnum, 2, p))
  varbetaAll = array(0, c(n_k, 2 * coefnum, 2 * coefnum, p))

  for (j in seq_along(k)) {
    for (i_p in 1:p) {
      betaCoefAll[j, , i_p] = beta(J, k[j], nu, Kernel_function, M_list[j], a = i_p, delta = 0)
      varbetaAll[j, , , i_p] = variance.estimator.v2(J, k[j], nu, Kernel_function, M_list[j], a = i_p, L, delta = 0)
      varbeta[j, , 1, i_p] = diag(varbetaAll[j, 1:coefnum, 1:coefnum, i_p])
      varbeta[j, , 2, i_p] = diag(varbetaAll[j, 1:coefnum + coefnum, 1:coefnum + coefnum, i_p])
    }
  }

  Test_tibble = NULL
  for (a in 1:p) {
    for (c in 1:p) {
      r_range = if (skeleton_only) 0 else 0:nu #(needs to be range adjusted later)
      for (r in r_range) {
        for (j in seq_along(k)) {
          if (!(a == c & r == 0)) {
            loc = r_to_loc(r, c, a, nu = nu, p = p)

            tmp1 = pnorm(abs(Re(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 1, a])), lower.tail = F) * 2
            tmp2 = pnorm(abs(Im(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 2, a])), lower.tail = F) * 2
            tmp3 = min(p.adjust(c(tmp1, tmp2), method = "BY"))

            tmp1_vec = c(Re(betaCoefAll[j, loc, a]), Im(betaCoefAll[j, loc, a]))
            matrix_tmp = varbetaAll[j, c(loc, loc + coefnum), c(loc, loc + coefnum), a]
            matrix_tmp[1, 2] = 0
            matrix_tmp[2, 1] = 0
            matrix_tmp = abs(matrix_tmp)
            tmp4 = pchisq(tmp1_vec %*% solve(matrix_tmp, b = tmp1_vec), lower.tail = F, df = 2)

            Test_tibble = tibble(
              Re = tmp1,
              Im = tmp2,
              padjust = tmp3,
              Chisq = tmp4,
              a = a,
              c = c,
              r = r,
              k = k[j],
              i = 1
            ) |>
              rbind(Test_tibble)
          }
        }
      }
    }
  }
  return(Test_tibble)
}

extract_skeleton = function(Edge_tibble, p, alpha = 0.05) {
  df = Edge_tibble |>
    filter(a < c, r == 0)

  edges = df |>
    group_by(a, c) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    ) |>
    filter(p_min < alpha)

  neighbors = vector("list", p)

  for (node in 1:p) {
    nb = unique(c(edges$c[edges$a == node], edges$a[edges$c == node]))
    neighbors[[node]] = sort(unique(c(node, nb)))
  }

  return(neighbors)
}

NonStGM_node_refit = function(x, neighbors, Kernel = 'Kernel_Triangular',
                              nu = 2, L = 1, alpha = 0.05) {

  Kernel_function = get(Kernel)

  Node_tibble = NULL

  for (a in seq_along(neighbors)) {
    S = neighbors[[a]]
    if (length(S) == 0) S = a
    x_sub = x[, S, drop = FALSE]
    a_local = which(S == a)
    # cat("dim(x_sub) =", dim(x_sub), "\n")
    # str(x_sub)
    variable_list = extractVariables(x_sub)
    J_sub = variable_list[[1]] ; p_sub = as.numeric(variable_list[2])
    coefnum_sub = (2 * p_sub * nu) + p_sub - 1
#
#     cat(
#       "\n=====================\n",
#       "Node =", a,
#       "\nS =", paste(S, collapse=","),
#       "\na_local =", a_local,
#       "\np_sub =", p_sub,
#       "\ncoefnum_sub =", coefnum_sub,
#       "\n=====================\n"
#     )

    k_sub = extractK(J_sub, coefnum_sub)

    n = nrow(x)
    M_grid = unique(round(seq(max(coefnum_sub , n^(1/5)),
                              max(coefnum_sub + 100 , n^(1/2)) , length.out = 10)))
    M_list = local_M_selection(J_sub, k_sub, M_grid)
    for (j in seq_along(k_sub)) {
      bc = beta(J_sub, k_sub[j], nu, Kernel_function, M_list[j], a = a_local, delta = 0)
      vb_full = variance.estimator.v2(J_sub, k_sub[j], nu, Kernel_function, M_list[j],
                                      a = a_local, L, delta = 0)
      vb_re = diag(vb_full[1:coefnum_sub, 1:coefnum_sub])
      vb_im = diag(vb_full[1:coefnum_sub + coefnum_sub, 1:coefnum_sub + coefnum_sub])

      for (r in 1:nu) {
        loc = r_to_loc(r, a_local, a_local, nu = nu, p = p_sub)

        tmp1 = pnorm(abs(Re(bc[loc])) / sqrt(abs(vb_re[loc])), lower.tail = F) * 2
        tmp2 = pnorm(abs(Im(bc[loc])) / sqrt(abs(vb_im[loc])), lower.tail = F) * 2
        tmp3 = min(p.adjust(c(tmp1, tmp2), method = "BY"))

        tmp1_vec = c(Re(bc[loc]), Im(bc[loc]))
        matrix_tmp = vb_full[c(loc, loc + coefnum_sub), c(loc, loc + coefnum_sub)]
        matrix_tmp[1, 2] = 0
        matrix_tmp[2, 1] = 0
        matrix_tmp = abs(matrix_tmp)
        tmp4 = pchisq(tmp1_vec %*% solve(matrix_tmp, b = tmp1_vec), lower.tail = F, df = 2)

        Node_tibble = tibble(
          Re = tmp1,
          Im = tmp2,
          padjust = tmp3,
          Chisq = tmp4,
          a = a,
          c = a,
          r = r,
          k = k_sub[j],
          i = 1
        ) |>
          rbind(Node_tibble)
      }
    }
  }
  return(Node_tibble)
}

graph_recovery_p5 = function(Test_tibble, alpha = 0.05) {
  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r != 0))
  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )
  true_edges = list(c(1,2), c(1,4), c(1,5), c(2,3), c(3,4), c(4,5))
  true_nodes = 1
  n_true_edges = length(true_edges)
  n_true_nodes = length(true_nodes)
  n_false_possible = (10 - n_true_edges) + (5 - n_true_nodes)
  recovery = Graph_tibble |>
    group_by(i) |>
    summarise(
      true_edges_found = sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      true_nodes_found = sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      false_edges = sum(a != c & p_min < alpha) - sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      false_nodes = sum(a == c & p_min < alpha) - sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      .groups = "drop"
    ) |>
    mutate(
      full = (true_edges_found == n_true_edges) & (true_nodes_found == n_true_nodes),
      edge_pct = true_edges_found / n_true_edges,
      node_pct = true_nodes_found / n_true_nodes,
      false_positives = pmax(0, false_edges) + pmax(0, false_nodes)
    )
  tibble(
    Method = "Global BY",
    Full_Detection = mean(recovery$full),
    Avg_Edge_Identification = mean(recovery$edge_pct),
    Avg_Node_Identification = mean(recovery$node_pct),
    False_Positive_Rate = if (n_false_possible > 0) mean(recovery$false_positives) / n_false_possible else 0
  )
}

graph_recovery_p5_twostep = function(Edge_tibble, Node_tibble, alpha = 0.05) {
  Edge_graph = Edge_tibble |>
    filter(a < c, r == 0) |>
    group_by(a, c, i) |>
    summarise(p_min = min(p.adjust(c(Re, Im), method = "BY")), .groups = "drop")
  Node_graph = Node_tibble |>
    group_by(a, c, i) |>
    summarise(p_min = min(p.adjust(c(Re, Im), method = "BY")), .groups = "drop")
  true_edges = list(c(1,2), c(1,4), c(1,5), c(2,3), c(3,4), c(4,5))
  true_nodes = 1
  n_true_edges = length(true_edges)
  n_true_nodes = length(true_nodes)
  n_false_possible = (10 - n_true_edges) + (5 - n_true_nodes)
  edge_recovery = Edge_graph |>
    group_by(i) |>
    summarise(
      true_edges_found = sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      false_edges = sum(p_min < alpha) - sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      .groups = "drop"
    )
  node_recovery = Node_graph |>
    group_by(i) |>
    summarise(
      true_nodes_found = sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      false_nodes = sum(a == c & p_min < alpha) - sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      .groups = "drop"
    )
  recovery = full_join(edge_recovery, node_recovery, by = "i") |>
    mutate(
      true_edges_found = coalesce(true_edges_found, 0L),
      false_edges = coalesce(false_edges, 0),
      true_nodes_found = coalesce(true_nodes_found, 0L),
      false_nodes = coalesce(false_nodes, 0),
      full = (true_edges_found == n_true_edges) & (true_nodes_found == n_true_nodes),
      edge_pct = true_edges_found / n_true_edges,
      node_pct = true_nodes_found / n_true_nodes,
      false_positives = pmax(0, false_edges) + pmax(0, false_nodes)
    )
  tibble(
    Method = "TwoStep BY",
    Full_Detection = mean(recovery$full),
    Avg_Edge_Identification = mean(recovery$edge_pct),
    Avg_Node_Identification = mean(recovery$node_pct),
    False_Positive_Rate = if (n_false_possible > 0) mean(recovery$false_positives) / n_false_possible else 0
  )
}

sim.tvVAR_p5_alltv = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,
    0.3,  0.6,  0,    0,    0,
    0,    0.3,  0.5,  0,    0,
    0,    0,    0.3, 0.55, 0,
    0.3, 0,    0,    0.3,  0.5
  ), ncol = 5, byrow = TRUE)
  n = m + burnin
  p = 5
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)
  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    A.t[2, 2] = st[tt]
    A.t[3, 3] = st[tt]
    A.t[4, 4] = st[tt]
    A.t[5, 5] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}

# graph_recovery_p5_alltv = function(Test_tibble, alpha = 0.05) {
#   Test_tibble = Test_tibble |>
#     filter(a <= c) |>
#     filter((a != c & r == 0) | (a == c & r != 0))
#
#   Graph_tibble = Test_tibble |>
#     group_by(a, c, i) |>
#     summarise(
#       p_min = min(p.adjust(c(Re, Im), method = "BY")),
#       .groups = "drop"
#     )
#
#   graph_recovery_results = Graph_tibble |>
#     group_by(i) |>
#     summarise(
#       detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
#       detect_22_tv = any(a == 2 & c == 2 & p_min < alpha),
#       detect_33_tv = any(a == 3 & c == 3 & p_min < alpha),
#       detect_44_tv = any(a == 4 & c == 4 & p_min < alpha),
#       detect_55_tv = any(a == 5 & c == 5 & p_min < alpha),
#       detect_12 = any(a == 1 & c == 2 & p_min < alpha),
#       detect_14 = any(a == 1 & c == 4 & p_min < alpha),
#       detect_15 = any(a == 1 & c == 5 & p_min < alpha),
#       detect_23 = any(a == 2 & c == 3 & p_min < alpha),
#       detect_34 = any(a == 3 & c == 4 & p_min < alpha),
#       detect_45 = any(a == 4 & c == 5 & p_min < alpha),
#       n_offdiag_edges = sum(a != c & p_min < alpha),
#       .groups = "drop"
#     )
#
#   accuracy_summary = graph_recovery_results |>
#     mutate(
#       all_offdiag_detected = detect_12 & detect_14 & detect_15 &
#         detect_23 & detect_34 & detect_45,
#       all_tv_detected = detect_11_tv & detect_22_tv & detect_33_tv &
#         detect_44_tv & detect_55_tv,
#       all_true_detected = all_offdiag_detected & all_tv_detected,
#       perfect = all_true_detected & (n_offdiag_edges == 6),
#       sensitivity = (detect_11_tv + detect_22_tv + detect_33_tv +
#                        detect_44_tv + detect_55_tv +
#                        detect_12 + detect_14 + detect_15 +
#                        detect_23 + detect_34 + detect_45) / 11,
#       average_tv_connection = (detect_11_tv + detect_22_tv + detect_33_tv +
#                                  detect_44_tv + detect_55_tv) / 5,
#       false_positives = pmax(0, n_offdiag_edges - 6)
#     )
#
#   tibble(
#     Method = "Global BY",
#     Perfect_Recovery = mean(accuracy_summary$perfect),
#     All_True_Detected = mean(accuracy_summary$all_true_detected),
#     Detect_11_TV = mean(accuracy_summary$detect_11_tv),
#     Detect_22_TV = mean(accuracy_summary$detect_22_tv),
#     Detect_33_TV = mean(accuracy_summary$detect_33_tv),
#     Detect_44_TV = mean(accuracy_summary$detect_44_tv),
#     Detect_55_TV = mean(accuracy_summary$detect_55_tv),
#     Average_TV_Connection = mean(accuracy_summary$average_tv_connection),
#     Detect_12 = mean(accuracy_summary$detect_12),
#     Detect_14 = mean(accuracy_summary$detect_14),
#     Detect_15 = mean(accuracy_summary$detect_15),
#     Detect_23 = mean(accuracy_summary$detect_23),
#     Detect_34 = mean(accuracy_summary$detect_34),
#     Detect_45 = mean(accuracy_summary$detect_45),
#     Mean_Sensitivity = mean(accuracy_summary$sensitivity),
#     Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
#     Mean_False_Positives = mean(accuracy_summary$false_positives)
#   )
# }
#
# graph_recovery_p5_alltv_twostep = function(Edge_tibble, Node_tibble, alpha = 0.05) {
#   Edge_graph = Edge_tibble |>
#     filter(a < c, r == 0) |>
#     group_by(a, c, i) |>
#     summarise(p_min = min(p.adjust(c(Re, Im), method = "BY")), .groups = "drop")
#
#   Node_graph = Node_tibble |>
#     group_by(a, c, i) |>
#     summarise(p_min = min(p.adjust(c(Re, Im), method = "BY")), .groups = "drop")
#
#   Combined = bind_rows(Edge_graph, Node_graph)
#
#   graph_recovery_results = Combined |>
#     group_by(i) |>
#     summarise(
#       detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
#       detect_22_tv = any(a == 2 & c == 2 & p_min < alpha),
#       detect_33_tv = any(a == 3 & c == 3 & p_min < alpha),
#       detect_44_tv = any(a == 4 & c == 4 & p_min < alpha),
#       detect_55_tv = any(a == 5 & c == 5 & p_min < alpha),
#       detect_12 = any(a == 1 & c == 2 & p_min < alpha),
#       detect_14 = any(a == 1 & c == 4 & p_min < alpha),
#       detect_15 = any(a == 1 & c == 5 & p_min < alpha),
#       detect_23 = any(a == 2 & c == 3 & p_min < alpha),
#       detect_34 = any(a == 3 & c == 4 & p_min < alpha),
#       detect_45 = any(a == 4 & c == 5 & p_min < alpha),
#       n_offdiag_edges = sum(a != c & p_min < alpha),
#       .groups = "drop"
#     )
#
#   accuracy_summary = graph_recovery_results |>
#     mutate(
#       all_offdiag_detected = detect_12 & detect_14 & detect_15 &
#         detect_23 & detect_34 & detect_45,
#       all_tv_detected = detect_11_tv & detect_22_tv & detect_33_tv &
#         detect_44_tv & detect_55_tv,
#       all_true_detected = all_offdiag_detected & all_tv_detected,
#       perfect = all_true_detected & (n_offdiag_edges == 6),
#       sensitivity = (detect_11_tv + detect_22_tv + detect_33_tv +
#                        detect_44_tv + detect_55_tv +
#                        detect_12 + detect_14 + detect_15 +
#                        detect_23 + detect_34 + detect_45) / 11,
#       average_tv_connection = (detect_11_tv + detect_22_tv + detect_33_tv +
#                                  detect_44_tv + detect_55_tv) / 5,
#       false_positives = pmax(0, n_offdiag_edges - 6)
#     )
#
#   tibble(
#     Method = "TwoStep BY",
#     Perfect_Recovery = mean(accuracy_summary$perfect),
#     All_True_Detected = mean(accuracy_summary$all_true_detected),
#     Detect_11_TV = mean(accuracy_summary$detect_11_tv),
#     Detect_22_TV = mean(accuracy_summary$detect_22_tv),
#     Detect_33_TV = mean(accuracy_summary$detect_33_tv),
#     Detect_44_TV = mean(accuracy_summary$detect_44_tv),
#     Detect_55_TV = mean(accuracy_summary$detect_55_tv),
#     Average_TV_Connection = mean(accuracy_summary$average_tv_connection),
#     Detect_12 = mean(accuracy_summary$detect_12),
#     Detect_14 = mean(accuracy_summary$detect_14),
#     Detect_15 = mean(accuracy_summary$detect_15),
#     Detect_23 = mean(accuracy_summary$detect_23),
#     Detect_34 = mean(accuracy_summary$detect_34),
#     Detect_45 = mean(accuracy_summary$detect_45),
#     Mean_Sensitivity = mean(accuracy_summary$sensitivity),
#     Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
#     Mean_False_Positives = mean(accuracy_summary$false_positives)
#   )
# }

graph_recovery_p5_alltv = function(Test_tibble, alpha = 0.05) {
  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r != 0))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  true_edges = list(c(1,2), c(1,4), c(1,5), c(2,3), c(3,4), c(4,5))
  true_nodes = 1:5
  n_true_edges = length(true_edges)
  n_true_nodes = length(true_nodes)
  n_false_possible = (10 - n_true_edges) + (5 - n_true_nodes)

  recovery = Graph_tibble |>
    group_by(i) |>
    summarise(
      true_edges_found = sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      true_nodes_found = sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      false_edges = sum(a != c & p_min < alpha) - sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      false_nodes = sum(a == c & p_min < alpha) - sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      .groups = "drop"
    ) |>
    mutate(
      full = (true_edges_found == n_true_edges) & (true_nodes_found == n_true_nodes),
      edge_pct = true_edges_found / n_true_edges,
      node_pct = true_nodes_found / n_true_nodes,
      false_positives = pmax(0, false_edges) + pmax(0, false_nodes)
    )

  tibble(
    Method = "Global BY",
    Full_Detection = mean(recovery$full),
    Avg_Edge_Identification = mean(recovery$edge_pct),
    Avg_Node_Identification = mean(recovery$node_pct),
    False_Positive_Rate = if (n_false_possible > 0) mean(recovery$false_positives) / n_false_possible else 0
  )
}

graph_recovery_p5_alltv_twostep = function(Edge_tibble, Node_tibble, alpha = 0.05) {
  Edge_graph = Edge_tibble |>
    filter(a < c, r == 0) |>
    group_by(a, c, i) |>
    summarise(p_min = min(p.adjust(c(Re, Im), method = "BY")), .groups = "drop")

  Node_graph = Node_tibble |>
    group_by(a, c, i) |>
    summarise(p_min = min(p.adjust(c(Re, Im), method = "BY")), .groups = "drop")

  true_edges = list(c(1,2), c(1,4), c(1,5), c(2,3), c(3,4), c(4,5))
  true_nodes = 1:5
  n_true_edges = length(true_edges)
  n_true_nodes = length(true_nodes)
  n_false_possible = (10 - n_true_edges) + (5 - n_true_nodes)

  edge_recovery = Edge_graph |>
    group_by(i) |>
    summarise(
      true_edges_found = sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      false_edges = sum(p_min < alpha) - sum(sapply(true_edges, function(e)
        any(a == e[1] & c == e[2] & p_min < alpha))),
      .groups = "drop"
    )

  node_recovery = Node_graph |>
    group_by(i) |>
    summarise(
      true_nodes_found = sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      false_nodes = sum(a == c & p_min < alpha) - sum(sapply(true_nodes, function(nd)
        any(a == nd & c == nd & p_min < alpha))),
      .groups = "drop"
    )

  recovery = full_join(edge_recovery, node_recovery, by = "i") |>
    mutate(
      true_edges_found = coalesce(true_edges_found, 0L),
      false_edges = coalesce(false_edges, 0),
      true_nodes_found = coalesce(true_nodes_found, 0L),
      false_nodes = coalesce(false_nodes, 0),
      full = (true_edges_found == n_true_edges) & (true_nodes_found == n_true_nodes),
      edge_pct = true_edges_found / n_true_edges,
      node_pct = true_nodes_found / n_true_nodes,
      false_positives = pmax(0, false_edges) + pmax(0, false_nodes)
    )

  tibble(
    Method = "TwoStep BY",
    Full_Detection = mean(recovery$full),
    Avg_Edge_Identification = mean(recovery$edge_pct),
    Avg_Node_Identification = mean(recovery$node_pct),
    False_Positive_Rate = if (n_false_possible > 0) mean(recovery$false_positives) / n_false_possible else 0
  )
}








######### V2
