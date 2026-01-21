sim.tvVAR = function(burnin, m, TV_size) {
  A = matrix(c(0.5, 0.2, 0, 0, 0.8, 0, 0, 0.3, 0.6), ncol = 3, byrow = T)
  n = m + burnin
  p = 3
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)
  for (tt in (2:n)) {
    A.t = A
    A.t[1, 1] = st[tt]  # Fixed indexing
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}
