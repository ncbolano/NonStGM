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
  k1 = (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 = sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }

  sigmaIm = diag(sigma1 + sigma2)

  return(sigmaIm)
}


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
  k1 = (k - n / 2 - 1)
  for (s1 in (-M:M)) {
    for (s2 in (-M:M)) {
      if (abs(s1 + s2 + 2 * k1) < L + 1) sigma2 = sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
    }
  }

  return(sigma1+sigma2)
}
