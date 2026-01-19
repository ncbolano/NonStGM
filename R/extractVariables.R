extractVariables = function(x) {
  p = ncol(x)
  nu = 2
  JJ0 = mvfft(x) / sqrt(nrow(x))
  n_rows = nrow(x)
  backward = c(((n_rows / 2) + 1):n_rows)
  JJ = rbind(JJ0[backward, ], JJ0[c(1:(n_rows / 2)), ])
  return(list(p,nu,JJ))
}

# Load from list
