#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
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
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
#'
whittle_loss = function(S_hat, I_k) {
  S_hat = Re(S_hat)

  # Regularization for invertibility if needed
  S_hat = S_hat + diag(1e-8, nrow(S_hat))

  return( sum(diag(solve(S_hat) %*% I_k)) + log(det(S_hat)) )
}
#' Transformation of frequencies (excluding a single index)
#'
#' @param J P-dimensional discrete fourier transform
#' @param k scalar value for \omega(k) = (2pi *k*) / n
#' @param nu scalar value which specifies amount of local DFT's to create frequency matrix *J_k^n*
#' @param a scalar index value which determines removal of specific index from local DFT matrix
#' @return
#' @noRd
#'
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
  S_k = matrix(0 + 0i, p, p)  # Complex matrix
  weight_sum = 0

  y = floor(n^(1/2))
  freq_window = (k - M):(k + M)

  if (k - M < 1) {
    JJ_smooth = rbind(JJ[((n-y):n),] , JJ) # Adding additional smoothing out of left edge case for the max M value of n^(1/2)
    for (j in freq_window) {
      weight = Kernel_func((k - j) / M)
      if (weight > 0) {
        I_j = outer(JJ_smooth[j+y, ], Conj(JJ_smooth[j+y, ]))
        S_k = S_k + weight * I_j
        weight_sum = weight_sum + weight
      }
    }
  }
  else {
    for (j in freq_window) {
      weight = Kernel_func((k - j) / M)
      if (weight > 0) {
        I_j = outer(JJ[j, ], Conj(JJ[j, ]))
        S_k = S_k + weight * I_j
        weight_sum = weight_sum + weight
      }
    }
  }

    # Normalize by sum of weights
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
  M_initial = 2 * coefnum
  n_freq = floor(nrow(JJ) / 2)

  trace_smooth = numeric(n_freq)
  diag_smooth = matrix(0, nrow = n_freq, ncol = ncol(JJ)) # Matrix for diagonal elements
  largest_eig_smooth = numeric(n_freq)
  condition_number_smooth = numeric(n_freq)

  for (k in 1:n_freq) {

    # Smoothed Periodogram (weighted by kernel)
    S_k = Re(smoothed_spectral_density(JJ, k, M_initial, Kernel_Triangular))
    trace_smooth[k] = sum(diag(S_k))
    diag_smooth[k, ] = diag(S_k)

    eigenvals_smooth = eigen(S_k, only.values = TRUE)$values
    eigenvals_smooth = Re(eigenvals_smooth)
    largest_eig_smooth[k] = max(eigenvals_smooth)
    condition_number_smooth[k] = max(eigenvals_smooth) / (min(eigenvals_smooth) + 1e-10)
  }

  composite_smooth = trace_smooth / condition_number_smooth # From frequency 1:(n/2)
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

      CV_scores[j] = whittle_loss(S_hat_k, I_candidate_k)

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
  M_grid = unique(round(seq(max(coefnum , n^(1/5)), n^(1/2), length.out = 50)))
  M_list = local_M_selection(JJ, k, M_grid)

  return(M_list)
}
# R = 500
# nu = 2
# p = 3
# Kernel = 'Kernel_Triangular'
# # alpha = .05
# burnin = 200
# n = 2^11
# TV_size = .6
#
# x = sim.tvVAR(burnin = 20, m = n, TV_size = TV_size)
# JJ0 = mvfft(x) / sqrt(nrow(x))
# n_rows = nrow(x)
# backward = c(((n_rows / 2) + 1):n_rows)
# JJ = rbind(JJ0[backward, ], JJ0[c(1:(n_rows / 2)), ])

# frequencies = (1:n_freq)
# plot_data_smooth = data.frame(
#   frequency = frequencies,
#   trace = Re(trace_smooth),
#   var1 = Re(diag_smooth[, 1]),
#   var2 = Re(diag_smooth[, 2]),
#   var3 = Re(diag_smooth[, 3]),
#   largest_eig = largest_eig_smooth,
#   cond_num = condition_number_smooth
# )
#
# p1_smooth = ggplot(plot_data_smooth, aes(x = frequency, y = trace)) +
#   geom_line(color = "darkgreen", linewidth = 0.8) +
#   labs(
#     title = "Trace of Smoothed Spectral Density Matrix",
#     subtitle = "Total power across all 3 variables (kernel-smoothed)",
#     x = "Frequency (cycles per observation)",
#     y = "Trace(S(ω))"
#   ) +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold"))
#
# print(p1_smooth)

# Need to append small amount at front and at end to account for max M testing size
# Need to ensure that M is at least c * coefnum

# M_initial = 2 * coefnum
#
# n_freq = floor(nrow(JJ) / 2)
#
# trace_smooth = numeric(n_freq)
# diag_smooth = matrix(0, nrow = n_freq, ncol = 3)
# largest_eig_smooth = numeric(n_freq)
# condition_number_smooth = numeric(n_freq)
#
# for (k in 1:n_freq) {
#
#   # Smoothed Periodogram (weighted by kernel)
#   S_k = smoothed_spectral_density(JJ, k, M_initial, Kernel_Triangular)
#   trace_smooth[k] = sum(diag(S_k))
#   diag_smooth[k, ] = diag(S_k)
#
#   eigenvals_smooth = eigen(S_k, only.values = TRUE)$values
#   eigenvals_smooth = Re(eigenvals_smooth)
#   largest_eig_smooth[k] = max(eigenvals_smooth)
#   condition_number_smooth[k] = max(eigenvals_smooth) / (min(eigenvals_smooth) + 1e-10)
# }
#
# composite_smooth = trace_smooth / condition_number_smooth # From frequency 1:(n/2)
# chosen_frequencies = return_max(composite_smooth , M_initial)
# k = chosen_frequencies
#
# # End of chosen list of two k frequencies
# coefnum = (2 * p * nu) + p - 1
# M_grid = unique(round(seq(max(coefnum , n^(1/5)), n^(1/2), length.out = 50)))
#
#
# CV_scores = numeric(length(M_grid))
# M_list = local_M_selection(JJ, k, M_grid)




