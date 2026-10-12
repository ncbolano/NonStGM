whittle_loss = function(S_hat, I_k) {
  S_hat = Re(S_hat)

  # Regularization for invertibility if needed
  S_hat = S_hat + diag(1e-8, nrow(S_hat))

  return( sum(diag(solve(S_hat) %*% I_k)) + log(det(S_hat)) )
}

smoothed_spectral_density_LOO = function(JJ, k, M, Kernel_func = Kernel_Triangular, leave_out = k) {
  n = nrow(JJ)
  p = ncol(JJ)
  S_k = matrix(0+0i, p, p)
  weight_sum = 0

  freq_window = max(1, k - M):min(n, k + M)

  for (j in freq_window) {
    if (j == leave_out) next
    weight = Kernel_func((k - j) / M)
    if (weight > 0) {
      I_j = outer(JJ[j,], Conj(JJ[j,]))
      S_k = S_k + weight * I_j
      weight_sum = weight_sum + weight
    }
  }

  if (weight_sum > 0) S_k = S_k / weight_sum
  return(S_k)
}

M_grid = unique(round(seq(2, 100, length.out = 25)))
n_freq = floor(nrow(JJ)/2)

CV_scores = numeric(length(M_grid))

for (m_index in seq_along(M_grid)) {
  M = M_grid[m_index]
  total_loss = 0

  for (k in 1:n_freq) {
    S_hat = smoothed_spectral_density_LOO(JJ, k, M, Kernel_Triangular, leave_out = ks[1])
    I_k   = outer(JJ[k,], Conj(JJ[k,]))

    total_loss = total_loss + whittle_loss(S_hat, I_k)
  }

  CV_scores[m_index] = total_loss
  cat("M =", M, " CV =", total_loss, "\n")
}

best_M = M_grid[ which.min(CV_scores) ]
cat("\nOptimal bandwidth:", best_M, "\n")

# Compute final estimator using best_M (no leave-out)
S_final = array(0, dim=c(p, p, n_freq))
for (k in 1:n_freq) {
  S_final[,,k] = smoothed_spectral_density(JJ, k, best_M, Kernel_Triangular)
}

M = best_M




local_M_selection = function(JJ, ks, M_grid) {

  best_local_M = numeric(length = length(ks))

  for (i in seq_along(ks)) {
    k_candidate = ks[i]

    I_candidate_k = outer(JJ[k_candidate, ], Conj(JJ[k_candidate, ]))

    CV_scores = numeric(length = length(M_grid))

    for (j in seq_along(M_grid)) {
      M = M_grid[j]

      S_hat_k = smoothed_spectral_density_LOO(JJ, k = k_candidate, M = M, leave_out = k_candidate)

      CV_scores[j] = whittle_loss(S_hat_k, I_candidate_k)

      cat("M =", M, "K =", k_candidate , " CV =", CV_scores[j], "\n")
    }

    # Find which M from the grid gave the minimum CV score for this k_candidate
    best_M_for_this_k = M_grid[which.min(CV_scores)]

    # Store the best M for this k
    best_local_M[i] = best_M_for_this_k
  }

  return(best_local_M)
}
M_list = local_M_selection(JJ, ks, M_grid)
M_list = c(30,40)
