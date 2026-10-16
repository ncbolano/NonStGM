# EEG_simulated.R -- self-contained p=12 EEG-like VAR simulation study.
# Source this entire file. You do not need full_function_list.R or other scripts.
# Only external dependency: tibble. Install once if needed: install.packages("tibble")
# stats, utils, graphics, and parallel are supplied with R; no library() calls needed.
#
# In interactive R/RStudio, source() runs the ONE-replication example at the bottom.
# To load functions only, first set:
# options(NonStGM.EEG_simulated.run_example = FALSE)
# To run the bottom example noninteractively, set that option to TRUE first.
# Change R=1L to R=100L in the bottom call for a repeated study, or call
# run_EEG_simulated(R=100L, n_cores=4L) after sourcing.
#
# Model: X[t] = A X[t-1] + D[t] z[t], z[t] iid N_12(0,I).
# A is three four-node blocks; spectral radius 0.4 and max absolute row sum 0.8.
# Innovation SD increases smoothly after sample 34000 of 50000, peaks at 4x,
# and settles to 2.5x baseline. The default mean is zero; estimation centers columns.
# True undirected edges: 1-2,2-3,2-4,5-6,6-7,6-8,9-10,10-11,10-12.
# Conditionally nonstationary nodes: 1,2,5,6,9,10.
# Time-varying edges: 1-2,5-6,9-10. Edge scores below measure skeleton recovery;
# they do not additionally classify each edge as stationary or time-varying.
# All twelve marginal variances change, despite six conditionally stationary nodes.
# TV_size=0 gives a stationary control with the same nine edges.
#
# Estimation: one-step only, simultaneous node and edge inference.
# Core copied from ncbolano/NonStGM R/full_function_list.R, blob
# 071e37980517aa4b4ef55de3cf3c3060bb529a5f (retrieved 2026-10-16).
# Numerical formulas and adaptive k/M selection are preserved. This copy adds
# periodic target indexing to match predictor indexing, failure checks, and k/M metadata.
# P5's BY adjustment is applied within EACH node/edge, not across the whole graph.
# Accuracy, sensitivity, specificity, false discoveries and exact recovery are reported.
# *_pct values are percentages, *_n are counts. With R=1, Monte Carlo SE is NA.
#
# All helper functions live in a private environment to avoid collisions with older
# versions already loaded in R. Public functions are exported below that block.

if (!requireNamespace("tibble", quietly = TRUE)) {
  stop('Package "tibble" is required. Run install.packages("tibble"), then source this file again.')
}

.EEG_simulated_module = local({
tibble = tibble::tibble

# Stable Gaussian VAR(1) with EEG-like time-varying innovation variance.
# A[row, column] is the coefficient from column at t-1 to row at t.

p12_eeg_spec = function() {
  B = matrix(c(
    0.40, 0,    0,    0,
    0.35, 0.40, 0,    0,
    0,    0.40, 0.40, 0,
    0,    0.40, 0,    0.40
  ), nrow = 4, byrow = TRUE)
  A = kronecker(diag(3), B)
  dimnames(A) = list(paste0("X", 1:12), paste0("X", 1:12))
  edges = data.frame(
    a = c(1L, 2L, 2L, 5L, 6L, 6L, 9L, 10L, 10L),
    c = c(2L, 3L, 4L, 6L, 7L, 8L, 10L, 11L, 12L),
    time_varying = rep(c(TRUE, FALSE, FALSE), 3)
  )
  ns_nodes = c(1L, 2L, 5L, 6L, 9L, 10L)
  adjacency = matrix(FALSE, 12, 12, dimnames = dimnames(A))
  adjacency[as.matrix(edges[c("a", "c")])] = TRUE
  adjacency = adjacency | t(adjacency)
  list(
    A = A, B = B, edges = edges, adjacency = adjacency,
    nonstationary_nodes = ns_nodes,
    stationary_nodes = setdiff(1:12, ns_nodes),
    nodes = data.frame(node = 1:12, conditionally_nonstationary = 1:12 %in% ns_nodes),
    spectral_radius = max(Mod(eigen(A, only.values = TRUE)$values)),
    max_abs_row_sum = max(rowSums(abs(A)))
  )
}

# C1 smooth ramp: exactly zero before 0 and exactly one after 1.
.p12_smoothstep = function(v) {
  v = pmax(0, pmin(1, v))
  v^2 * (3 - 2 * v)
}

# A plain numeric m-by-12 matrix is returned by default, like sim.tvVAR_p15().
# TV_size is the ADDITIONAL peak innovation SD multiplier, not an AR coefficient:
# TV_size=3 gives baseline 1, peak 4, and final multiplier 2.5 by default.
# The onset is measured in retained samples, AFTER removing burnin.
# mean_level=0 is recommended for estimation; mean_level=-250 gives an EEG-like offset.
# All twelve marginal series change variance; only six nodes are conditionally NS.
sim.tvVAR_p12_eeg = function(burnin = 1000L, m = 50000L, TV_size = 3,
                            sigma0 = 25, mean_level = 0, seed = NULL,
                            onset = 0.68, rise_width = 0.08,
                            fall_start = 0.88, fall_width = 0.06,
                            late_fraction = 0.5, return_details = FALSE) {
  scalar_finite = function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar_finite(m) || m < 20 || m != floor(m)) stop("m must be an integer >= 20.")
  if (!scalar_finite(burnin) || burnin < 0 || burnin != floor(burnin))
    stop("burnin must be a nonnegative integer.")
  if (!scalar_finite(TV_size) || TV_size < 0) stop("TV_size must be nonnegative.")
  pars = c(onset, rise_width, fall_start, fall_width, late_fraction)
  if (length(pars) != 5L || any(!is.finite(pars)) ||
      onset <= 0 || rise_width <= 0 || fall_width <= 0 ||
      onset + rise_width > fall_start || fall_start + fall_width > 1 ||
      late_fraction < 0 || late_fraction > 1)
    stop("Require 0 < onset < onset+rise_width <= fall_start < fall_start+fall_width <= 1; late_fraction in [0,1].")
  if (!is.numeric(sigma0) || !length(sigma0) %in% c(1L, 12L) ||
      any(!is.finite(sigma0)) || any(sigma0 <= 0))
    stop("sigma0 must contain one or twelve positive finite innovation SDs.")
  if (!is.numeric(mean_level) || !length(mean_level) %in% c(1L, 12L) ||
      any(!is.finite(mean_level))) stop("mean_level must contain one or twelve finite means.")
  if (!is.logical(return_details) || length(return_details) != 1L || is.na(return_details))
    stop("return_details must be TRUE or FALSE.")
  if (!is.null(seed)) set.seed(seed)
  sigma0 = rep(sigma0, length.out = 12L)
  mu = rep(mean_level, length.out = 12L)
  spec = p12_eeg_spec()
  A = spec$A
  stopifnot(spec$spectral_radius < 1, spec$max_abs_row_sum < 1)

  # Exact stationary initialization in the baseline regime: P = A P A' + Sigma0.
  # Subsequent burn-in also stays in this baseline regime.
  P0 = matrix(solve(diag(144) - kronecker(A, A),
                    as.vector(diag(sigma0^2))), 12, 12)
  P0 = (P0 + t(P0)) / 2
  state = drop(t(chol(P0)) %*% rnorm(12))

  n = as.integer(m + burnin)
  retained_time = seq_len(n) - burnin
  u = retained_time / m
  h = .p12_smoothstep((u - onset) / rise_width) *
    (1 - (1 - late_fraction) * .p12_smoothstep((u - fall_start) / fall_width))
  g = 1 + TV_size * h
  z = matrix(rnorm(n * 12), nrow = n, ncol = 12)
  y = matrix(0, nrow = n, ncol = 12)
  for (tt in seq_len(n)) {
    sd_t = sigma0
    sd_t[spec$nonstationary_nodes] = sd_t[spec$nonstationary_nodes] * g[tt]
    state = drop(A %*% state) + sd_t * z[tt, ]
    y[tt, ] = state
  }
  keep = burnin + seq_len(m)
  x = sweep(y[keep, , drop = FALSE], 2, mu, "+")
  colnames(x) = paste0("X", 1:12)
  if (!return_details) return(x)

  # TV_size=0 is a stationary null-control variant, with the same edge skeleton.
  if (TV_size == 0) {
    spec$nonstationary_nodes = integer(0)
    spec$stationary_nodes = 1:12
    spec$nodes$conditionally_nonstationary = FALSE
    spec$edges$time_varying = FALSE
  }
  list(
    x = x, A = A, truth = spec,
    time = seq_len(m), sd_multiplier = g[keep],
    settings = list(m = m, burnin = burnin, TV_size = TV_size,
                    sigma0 = sigma0, mean_level = mu, seed = seed,
                    onset = onset, rise_width = rise_width, fall_start = fall_start,
                    fall_width = fall_width, late_fraction = late_fraction),
    transition_samples = c(rise_start = m * onset, rise_end = m * (onset + rise_width),
                           fall_start = m * fall_start, fall_end = m * (fall_start + fall_width))
  )
}

# Exact interior precision blocks for the Gaussian VAR trajectory.
# This is not the finite-sample marginal inverse covariance at its endpoints.
p12_precision_blocks = function(sample_index, m = 50000L, TV_size = 3, sigma0 = 25,
                               onset = 0.68, rise_width = 0.08,
                               fall_start = 0.88, fall_width = 0.06,
                               late_fraction = 0.5) {
  spec = p12_eeg_spec()
  sigma0 = rep(sigma0, length.out = 12)
  Q = function(time) {
    u = time / m
    h = .p12_smoothstep((u - onset) / rise_width) *
      (1 - (1 - late_fraction) * .p12_smoothstep((u - fall_start) / fall_width))
    s = sigma0
    s[spec$nonstationary_nodes] = s[spec$nonstationary_nodes] * (1 + TV_size * h)
    diag(1 / s^2)
  }
  list(same_time = Q(sample_index) + t(spec$A) %*% Q(sample_index + 1) %*% spec$A,
       lag_one = -Q(sample_index) %*% spec$A)
}

# One-step estimation core.
extractVariables = function(x) {
    p = ncol(x)
    nu = 2
    JJ0 = mvfft(x)/sqrt(nrow(x))
    n_rows = nrow(x)
    backward = c(((n_rows/2) + 1):n_rows)
    JJ = rbind(JJ0[backward, , drop = FALSE], JJ0[1:(n_rows/2), , drop = FALSE])
    return(list(JJ, p, nu))
}

J.to.J.k.nu.a = function(J, k, nu, a) {
    n = nrow(J)
    p = ncol(J)
    indices_to_get = k + (-nu:nu)
    wrapped_indices = ((indices_to_get - 1)%%n) + 1
    Jtemp = c(t(J[wrapped_indices, ]))
    Jtemp2 = Jtemp[-(p * nu + a)]
    return(Jtemp2)
}

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

whittle_loss = function(S_hat, I_k) {
    S_hat = Re(S_hat)
    S_hat = S_hat + diag(1e-06, nrow(S_hat))
    return(Re(sum(diag(solve(S_hat) %*% I_k)) + log(det(S_hat))))
}

smoothed_spectral_density_LOO = function(JJ, k, M, Kernel_func = Kernel_Triangular, leave_out = k) {
    n = nrow(JJ)
    p = ncol(JJ)
    S_k = matrix(0 + (0+0i), p, p)
    weight_sum = 0
    for (j_offset in -M:M) {
        current_j = k + j_offset
        idx = ((current_j - 1)%%n) + 1
        if (idx == leave_out) 
            next
        weight = Kernel_func(j_offset/M)
        if (weight > 0) {
            I_j = outer(JJ[idx, ], Conj(JJ[idx, ]))
            S_k = S_k + weight * I_j
            weight_sum = weight_sum + weight
        }
    }
    if (weight_sum > 0) 
        S_k = S_k/weight_sum
    return(S_k)
}

smoothed_spectral_density = function(JJ, k, M, Kernel_func = Kernel_Triangular) {
    n = nrow(JJ)
    p = ncol(JJ)
    S_k = matrix(0 + (0+0i), p, p)
    weight_sum = 0
    for (j_offset in -M:M) {
        current_j = k + j_offset
        idx = ((current_j - 1)%%n) + 1
        weight = Kernel_func(j_offset/M)
        if (weight > 0) {
            I_j = outer(JJ[idx, ], Conj(JJ[idx, ]))
            S_k = S_k + weight * I_j
            weight_sum = weight_sum + weight
        }
    }
    if (weight_sum > 0) {
        S_k = S_k/weight_sum
    }
    return(S_k)
}

return_max = function(x, M) {
    max1 = which.max(x)
    x2 = x
    min_dist = M
    forbidden_indices = max(1, max1 - min_dist):min(length(x), max1 + min_dist)
    x2[forbidden_indices] = -Inf
    max2 = which.max(x2)
    if (is.infinite(x2[max2])) {
        warning("No valid second maximum found at distance >= 2*M")
        return(c(max1, NA))
    }
    return(c(max1, max2))
}

extractK = function(JJ, coefnum) {
    cat("\n--- extractK ---\n")
    cat("dim(JJ):", paste(dim(JJ), collapse = " x "), "\n")
    cat("coefnum:", coefnum, "\n")
    M_initial = 2 * coefnum
    n_freq = floor(nrow(JJ)/2)
    floor = M_initial + 1
    ceiling = n_freq - M_initial
    cat("M_initial =", M_initial, "\n")
    cat("n_freq =", n_freq, "\n")
    cat("floor =", floor, "\n")
    cat("ceiling =", ceiling, "\n")
    trace_smooth = numeric(n_freq)
    diag_smooth = matrix(0, nrow = n_freq, ncol = ncol(JJ))
    largest_eig_smooth = numeric(n_freq)
    condition_number_smooth = numeric(n_freq)
    for (k in floor:ceiling) {
        S_k = Re(smoothed_spectral_density(JJ, k, M_initial, Kernel_Triangular))
        trace_smooth[k] = sum(diag(S_k))
        diag_smooth[k, ] = diag(S_k)
        eigenvals_smooth = eigen(S_k, only.values = TRUE)$values
        eigenvals_smooth = Re(eigenvals_smooth)
        largest_eig_smooth[k] = max(eigenvals_smooth)
        condition_number_smooth[k] = max(eigenvals_smooth)/(min(eigenvals_smooth) + 1e-10)
    }
    composite_smooth = trace_smooth/((condition_number_smooth)^(0.5))
    chosen_frequencies = return_max(composite_smooth, M_initial)
    k = chosen_frequencies
    return(k)
}

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
        }
        best_M_for_this_k = M_grid[which.min(CV_scores)]
        best_local_M[i] = best_M_for_this_k
    }
    print(best_local_M)
    return(best_local_M)
}

Kernel_Uniform = function(x) {
    return(as.numeric(abs(x) <= 1))
}

Hat.R.reduced = function(J, k, nu, W, M1, a) {
    n = nrow(J)
    p = ncol(J)
    R_nu_M = matrix(0, (p * (2 * nu + 1) - 1), (p * (2 * nu + 1) - 1))
    for (s in (-M1:M1)) {
        J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
        R_nu_M = R_nu_M + W(s/M1) * J_k_S %*% t(Conj(J_k_S))
    }
    norm = sum(sapply((-M1:M1), function(s) W(s/M1)))
    return(R_nu_M/norm)
}

Hat.r = function(J, k, nu, W, M, a) {
    n = nrow(J)
    p = ncol(J)
    r_nu_M = rep(0, (p * (2 * nu + 1) - 1))
    for (s in (-M:M)) {
        J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
        r_nu_M = r_nu_M + W(s/M) * J_k_S * Conj(J[((k + s - 1) %% nrow(J)) + 1, a])
    }
    norm = sum(sapply((-M:M), function(s) W(s/M)))
    return(r_nu_M/norm)
}

Kernel_Triangular = function(x) {
    return(as.numeric(abs(x) <= 1) * (1 - abs(x)))
}

Kernel_Quadratic = function(x) {
    return(as.numeric(abs(x) <= 1) * (1 - abs(x)^2))
}

beta = function(J, k, nu, W, M, a, delta) {
    HatR = Hat.R.reduced(J, k, nu, W, M, a)
    if (delta > 0) {
        dim1 = ncol(HatR)
        HatR = HatR + diag(rep(delta, dim1))
    }
    if (rcond(HatR) < .Machine$double.eps) {
        p = ncol(J)
        len_out = p * (2 * nu + 1) - 1
        return(rep(NA, len_out))
    }
    hatr = Hat.r(J, k, nu, W, M, a)
    betaC = solve(a = HatR, b = hatr)
    return(beta_beta_r(betaC, nu, ncol(J)))
}

beta_beta_r = function(beta, nu, p) {
    tmp1 = as.vector(unlist(sapply((-nu):nu, function(i) {
        if (i == 0) 
            return((nu) * p + seq_len(p - 1))
        if (i < 0) 
            return((abs(i) + nu) * p + 1:p - 1)
        if (i > 0) 
            return((nu - abs(i)) * p + 1:p)
    })))
    tmp2 = seq_along(beta)
    beta_erg = beta
    beta_erg[tmp2 != tmp1] = ((beta + Conj(beta[tmp1]))[tmp2 != tmp1])/2
    return(beta_erg)
}

variance.estimator.v2 = function(J, k, nu, W, M, a, L, delta) {
    betaC = beta(J, k, nu, W, M, a, delta)
    if (any(!is.finite(betaC))) stop(sprintf("Singular regression at node %d, frequency %d, bandwidth %d.", a, k, M))
    dim1 = length(betaC)
    HatRred = Hat.R.reduced(J, k, nu, W, M, a)
    ridge = diag(rep(delta, dim1))
    invHatRreduced = solve(HatRred + ridge)
    norm = sum(sapply((-M:M), function(s) W(s/M)))
    Y = matrix(rep(0, 2 * (2 * M + 1) * dim1), nrow = 2 * dim1)
    p = ncol(J)
    for (s in (-M:M)) {
        J_k_S = J.to.J.k.nu.a(J, k + s, nu, a)
        temp1 = (Conj(J[((k + s - 1) %% nrow(J)) + 1, a]) - sum(betaC * Conj(J_k_S)))
        temp2 = invHatRreduced %*% J_k_S
        Y[, s + M + 1] = W(s/M) * c(Re(temp1 * temp2/norm), Im(temp1 * temp2/norm))
    }
    sigma1 = matrix(rep(0, (2 * dim1)^2), ncol = 2 * dim1)
    for (s1 in (1:(2 * M + 1))) {
        for (s2 in (1:(2 * M + 1))) {
            if (abs(s1 - s2) < L + 1) 
                sigma1 = sigma1 + Y[, s1] %*% t(Y[, s2])
        }
    }
    sigma2 = matrix(rep(0, (2 * dim1)^2), ncol = 2 * dim1)
    n = nrow(J)
    k1 = (k - n/2 - 1)
    for (s1 in (-M:M)) {
        for (s2 in (-M:M)) {
            if (abs(s1 + s2 + 2 * k1) < L + 1) 
                sigma2 = sigma2 + Y[, s1 + M + 1] %*% t(Y[, s2 + M + 1])
        }
    }
    return(sigma1 + sigma2)
}

NonStGM_adj = function(x, Kernel = "Kernel_Triangular", nu = 2, L = 1, alpha = 0.05, skeleton_only = FALSE) {
    Kernel_function = get(Kernel)
    variable_list = extractVariables(x)
    J = variable_list[[1]]
    p = as.numeric(variable_list[2])
    coefnum = (2 * p * nu) + p - 1
    k = extractK(J, coefnum)
    if (length(k) != 2L || any(!is.finite(k)) || anyDuplicated(k))
        stop("Frequency selection did not return two distinct finite bins.")
    n = nrow(x)
    M_grid = unique(round(seq(max(coefnum, n^(1/5)), max(coefnum + 100, n^(1/2)), length.out = 10)))
    M_list = local_M_selection(J, k, M_grid)
    n_k = length(k)
    betaCoefAll = array(0, c(n_k, coefnum, p))
    varbeta = array(0, c(n_k, coefnum, 2, p))
    varbetaAll = array(0, c(n_k, 2 * coefnum, 2 * coefnum, p))
    for (j in seq_along(k)) {
        for (i_p in 1:p) {
            betaCoefAll[j, , i_p] = beta(J, k[j], nu, Kernel_function, M_list[j], a = i_p, delta = 0)
            varbetaAll[j, , , i_p] = variance.estimator.v2(J, k[j], nu, Kernel_function, M_list[j], a = i_p, L, 
                delta = 0)
            varbeta[j, , 1, i_p] = diag(varbetaAll[j, 1:coefnum, 1:coefnum, i_p])
            varbeta[j, , 2, i_p] = diag(varbetaAll[j, 1:coefnum + coefnum, 1:coefnum + coefnum, i_p])
        }
    }
    Test_tibble = NULL
    for (a in 1:p) {
        for (c in 1:p) {
            r_range = if (skeleton_only) 
                0
            else 0:nu
            for (r in r_range) {
                for (j in seq_along(k)) {
                  if (!(a == c & r == 0)) {
                    loc = r_to_loc(r, c, a, nu = nu, p = p)
                    tmp1 = pnorm(abs(Re(betaCoefAll[j, loc, a]))/sqrt(abs(varbeta[j, loc, 1, a])), lower.tail = F) * 
                      2
                    tmp2 = pnorm(abs(Im(betaCoefAll[j, loc, a]))/sqrt(abs(varbeta[j, loc, 2, a])), lower.tail = F) * 
                      2
                    tmp3 = min(p.adjust(c(tmp1, tmp2), method = "BY"))
                    tmp1_vec = c(Re(betaCoefAll[j, loc, a]), Im(betaCoefAll[j, loc, a]))
                    matrix_tmp = varbetaAll[j, c(loc, loc + coefnum), c(loc, loc + coefnum), a]
                    matrix_tmp[1, 2] = 0
                    matrix_tmp[2, 1] = 0
                    matrix_tmp = abs(matrix_tmp)
                    tmp4 = pchisq(tmp1_vec %*% solve(matrix_tmp, b = tmp1_vec), lower.tail = F, df = 2)
                    Test_tibble = rbind(tibble(Re = tmp1, Im = tmp2, padjust = tmp3, Chisq = tmp4, a = a, c = c, 
                      r = r, k = k[j], i = 1), Test_tibble)
                  }
                }
            }
        }
    }
    attr(Test_tibble, "selected_frequencies") = k
    attr(Test_tibble, "selected_bandwidths") = M_list
    return(Test_tibble)
}


# One-step scoring follows the latest graph_recovery_p5():
# edges use a < c, r == 0; nodes use a == c, r > 0.
# BY is applied separately within EACH edge/node across its Re/Im tests.
# This reproduces the existing decision rule; it is not a global graph FDR correction.
# Rates below are percentages on [0,100]; counts end in _n.

.eeg_binary_metrics = function(truth, detected, prefix) {
  tp = sum(truth & detected)
  tn = sum(!truth & !detected)
  fp = sum(!truth & detected)
  fn = sum(truth & !detected)
  pct = function(num, den) if (den == 0) NA_real_ else 100 * num / den
  sensitivity = pct(tp, tp + fn)
  specificity = pct(tn, tn + fp)
  values = c(tp_n = tp, tn_n = tn, fp_n = fp, fn_n = fn,
             accuracy_pct = pct(tp + tn, length(truth)),
             sensitivity_pct = sensitivity, specificity_pct = specificity,
             precision_pct = pct(tp, tp + fp),
             fpr_pct = pct(fp, fp + tn),
             # The realized false-discovery proportion is zero with no discoveries.
             # Its average across replications estimates the false-discovery rate.
             fdr_pct = if (tp + fp == 0) 0 else pct(fp, tp + fp),
             balanced_accuracy_pct = (sensitivity + specificity) / 2)
  stats::setNames(values, paste0(prefix, "_", names(values)))
}

# Score ONE replication only. Pass nu and the estimator's selected frequencies
# to detect an omitted harmonic/frequency, as well as omitted node/edge decisions.
score_EEG_p12 = function(tests, truth, alpha = 0.05, nu = NULL,
                         frequencies = NULL) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must be strictly between zero and one.")
  required = c("a", "c", "r", "k", "Re", "Im")
  if (!is.data.frame(tests) || !all(required %in% names(tests)) || !nrow(tests))
    stop("tests must be a nonempty data frame with a, c, r, k, Re, and Im columns.")
  if (!all(vapply(tests[required], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(tests[required]))))
    stop("Test identifiers and p-values must be finite numeric values.")
  if (any(tests$Re < 0 | tests$Re > 1 | tests$Im < 0 | tests$Im > 1))
    stop("Every Re and Im p-value must lie in [0,1].")
  ids = as.matrix(tests[c("a", "c", "r", "k")])
  if (any(ids != floor(ids)) || any(tests$a < 1 | tests$a > 12) ||
      any(tests$c < 1 | tests$c > 12) || any(tests$r < 0) || any(tests$k < 0))
    stop("Invalid node, harmonic, or frequency identifiers.")
  if ("i" %in% names(tests) && (anyNA(tests$i) || length(unique(tests$i)) != 1L))
    stop("score_EEG_p12() accepts one replication at a time.")
  if (anyDuplicated(tests[c("a", "c", "r", "k")]))
    stop("Duplicate ordered (a,c,r,k) tests; repetitions must not be pooled.")
  if (is.null(nu)) nu = max(tests$r[tests$a == tests$c], 0)
  if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) ||
      nu < 1 || nu != floor(nu) || any(tests$r > nu))
    stop("nu must be a positive integer covering all tested harmonics.")
  if (is.null(frequencies)) frequencies = sort(unique(tests$k))
  if (!is.numeric(frequencies) || !length(frequencies) ||
      any(!is.finite(frequencies)) || any(frequencies != floor(frequencies)) ||
      any(frequencies < 0) || anyDuplicated(frequencies) ||
      !setequal(tests$k, frequencies))
    stop("frequencies must uniquely identify all selected frequency bins.")

  if (!is.list(truth) || !is.data.frame(truth$edges) ||
      !all(c("a", "c") %in% names(truth$edges)) ||
      !is.data.frame(truth$nodes) ||
      !all(c("node", "conditionally_nonstationary") %in% names(truth$nodes)))
    stop("truth must contain edges and nodes tables from sim.tvVAR_p12_eeg().")
  te = truth$edges[c("a", "c")]
  tn = truth$nodes
  if (!all(vapply(te, is.numeric, logical(1))) || any(!is.finite(as.matrix(te))) ||
      any(as.matrix(te) != floor(as.matrix(te))) ||
      any(te$a < 1 | te$c > 12 | te$a >= te$c) || anyDuplicated(te))
    stop("Truth edges must be unique unordered pairs with 1 <= a < c <= 12.")
  if (!is.numeric(tn$node) || !setequal(tn$node, 1:12) ||
      nrow(tn) != 12L || anyNA(tn$node) || anyDuplicated(tn$node) ||
      !is.logical(tn$conditionally_nonstationary) ||
      anyNA(tn$conditionally_nonstationary))
    stop("Truth must provide exactly one logical nonstationarity label per node 1:12.")

  pairs = as.data.frame(t(utils::combn(1:12, 2)))
  names(pairs) = c("a", "c")
  relevant = tests[(tests$a < tests$c & tests$r == 0) |
                     (tests$a == tests$c & tests$r > 0), , drop = FALSE]
  key = function(a, c, r, k) paste(a, c, r, k, sep = ":")
  edge_expected = expand.grid(pair = seq_len(nrow(pairs)), k = frequencies)
  node_expected = expand.grid(a = 1:12, r = seq_len(nu), k = frequencies)
  expected = c(key(pairs$a[edge_expected$pair], pairs$c[edge_expected$pair],
                   0, edge_expected$k),
               key(node_expected$a, node_expected$a, node_expected$r, node_expected$k))
  observed = key(relevant$a, relevant$c, relevant$r, relevant$k)
  if (length(observed) != length(expected) || !setequal(observed, expected))
    stop("Incomplete edge/node tests: require every edge at r=0 and every node at r=1:nu at every selected frequency.")

  decision_p = function(a, c) {
    z = relevant[relevant$a == a & relevant$c == c, , drop = FALSE]
    min(stats::p.adjust(c(z$Re, z$Im), method = "BY"))
  }
  edges = pairs
  edges$truth = paste(edges$a, edges$c) %in% paste(te$a, te$c)
  edges$p_min_BY = mapply(decision_p, edges$a, edges$c)
  edges$detected = edges$p_min_BY < alpha
  edges$correct = edges$detected == edges$truth
  nodes = data.frame(node = 1:12,
                     truth = tn$conditionally_nonstationary[match(1:12, tn$node)])
  nodes$p_min_BY = vapply(1:12, function(a) decision_p(a, a), numeric(1))
  nodes$detected = nodes$p_min_BY < alpha
  nodes$correct = nodes$detected == nodes$truth
  all_truth = c(edges$truth, nodes$truth)
  all_detected = c(edges$detected, nodes$detected)
  metrics = c(.eeg_binary_metrics(edges$truth, edges$detected, "edge"),
              .eeg_binary_metrics(nodes$truth, nodes$detected, "node"),
              .eeg_binary_metrics(all_truth, all_detected, "combined"),
              edge_graph_exact_pct = 100 * all(edges$correct),
              node_labels_exact_pct = 100 * all(nodes$correct),
              full_graph_exact_pct = 100 * all(c(edges$correct, nodes$correct)),
              # Full detection allows false positives; exact graph recovery does not.
              full_detection_pct = 100 * all(all_detected[all_truth]))
  list(edge_decisions = edges, node_decisions = nodes,
       metrics = as.data.frame(as.list(metrics)),
       alpha = alpha, nu = nu, frequencies = frequencies)
}

# Summaries use all completed replications. Undefined sensitivity/specificity or
# precision remains NA; n_valid shows the denominator for each metric's average.
# A single replication has NA Monte Carlo SE (it cannot estimate simulation error).
summarize_EEG_p12 = function(scores) {
  if (!is.list(scores) || !length(scores) ||
      !all(vapply(scores, function(z) is.list(z) &&
                    all(c("edge_decisions", "node_decisions", "metrics") %in% names(z)),
                  logical(1)))) stop("scores must be a nonempty list of score_EEG_p12() results.")
  reference = scores[[1]]
  for (z in scores) {
    if (!identical(z$edge_decisions[c("a", "c", "truth")],
                   reference$edge_decisions[c("a", "c", "truth")]) ||
        !identical(z$node_decisions[c("node", "truth")],
                   reference$node_decisions[c("node", "truth")]) ||
        !identical(names(z$metrics), names(reference$metrics)) ||
        !identical(z$alpha, reference$alpha) ||
        !identical(z$nu, reference$nu))
      stop("Replication summaries require the same truth, decisions, metrics, and alpha.")
  }
  per_rep = do.call(rbind, lapply(scores, `[[`, "metrics"))
  rownames(per_rep) = NULL
  summary_rows = lapply(names(per_rep), function(metric) {
    x = per_rep[[metric]]
    x = x[is.finite(x)]
    n = length(x)
    data.frame(metric = metric,
               unit = if (endsWith(metric, "_pct")) "percent" else "count",
               mean = if (n) mean(x) else NA_real_,
               mcse = if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_,
               n_valid = n)
  })
  frequencies = function(field, id_cols) {
    result = reference[[field]][c(id_cols, "truth")]
    detected = do.call(cbind, lapply(scores, function(z) z[[field]]$detected))
    correct = do.call(cbind, lapply(scores, function(z) z[[field]]$correct))
    result$n_replications = ncol(detected)
    result$detected_n = rowSums(detected)
    result$detection_pct = 100 * rowMeans(detected)
    result$accuracy_pct = 100 * rowMeans(correct)
    result
  }
  list(per_replication_metrics = cbind(replication = seq_along(scores), per_rep),
       metric_summary = do.call(rbind, summary_rows),
       edge_detection_frequency = frequencies("edge_decisions", c("a", "c")),
       node_detection_frequency = frequencies("node_decisions", "node"))
}

# Public fitting wrapper for an existing p=12 data matrix.
# Nodes and edges are estimated together with skeleton_only=FALSE.
fit_EEG_simulated = function(x, nu = 2L, L = 1L,
                             Kernel = "Kernel_Triangular", verbose = TRUE) {
  .eeg_check_integer(nu, "nu", 1L)
  .eeg_check_integer(L, "L", 0L)
  .eeg_check_flag(verbose, "verbose")
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) != 12L ||
      any(!is.finite(x))) stop("x must be a finite numeric matrix with 12 columns.")
  coefnum = 12 * (2 * nu + 1) - 1
  # Guarantees two eligible frequency bins whatever the first selected bin is.
  min_n = 16 * coefnum + 4
  if (nrow(x) %% 2 != 0 || nrow(x) < min_n)
    stop(sprintf("Use an even observation count of at least %d for nu=%d.", min_n, nu))
  if (!is.character(Kernel) || length(Kernel) != 1L || is.na(Kernel) ||
      !Kernel %in% c("Kernel_Triangular", "Kernel_Uniform", "Kernel_Quadratic"))
    stop("Kernel must name Kernel_Triangular, Kernel_Uniform, or Kernel_Quadratic.")
  if (any(apply(x, 2, stats::sd) == 0)) stop("Each column of x must vary.")
  x = sweep(x, 2, colMeans(x), "-")
  # Capturing diagnostic prints does not suppress warnings or errors.
  diagnostic_log = utils::capture.output(
    tests <- NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L,
                        skeleton_only = FALSE)
  )
  if (any(!is.finite(as.matrix(tests[c("Re", "Im")]))))
    stop("The one-step estimator produced nonfinite p-values; inspect the regression fit.")
  attr(tests, "diagnostic_log") = diagnostic_log
  if (verbose) cat(paste(diagnostic_log, collapse = "\n"), "\n")
  tests
}

.eeg_check_integer = function(x, name, lower = 0L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > .Machine$integer.max)
    stop(sprintf("%s must be an integer in [%d, %d].", name, lower, .Machine$integer.max))
}
.eeg_check_flag = function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(paste(name, "must be TRUE or FALSE."))
}

# All arguments needed by a worker are explicit so this also works on Windows PSOCK.
.eeg_run_one = function(replication, simulation_args, fit_args, seeds, alpha,
                       keep_tests, verbose) {
  tryCatch({
    if (verbose) message(sprintf("Replication %d (seed %d): simulating and fitting...",
                                 replication, seeds[replication]))
    start = proc.time()[["elapsed"]]
    simulation_args$seed = seeds[replication]
    sim = do.call(sim.tvVAR_p12_eeg, simulation_args)
    tests = do.call(fit_EEG_simulated, c(list(x = sim$x, verbose = FALSE), fit_args))
    tests$i = replication
    score = score_EEG_p12(tests, sim$truth, alpha = alpha, nu = fit_args$nu,
                          frequencies = attr(tests, "selected_frequencies"))
    elapsed = proc.time()[["elapsed"]] - start
    if (verbose) message(sprintf("Replication %d complete in %.1f seconds.", replication, elapsed))
    list(score = score, seed = seeds[replication], elapsed_seconds = elapsed,
         frequencies = attr(tests, "selected_frequencies"),
         bandwidths = attr(tests, "selected_bandwidths"),
         simulation = if (replication == 1L) sim else NULL,
         tests = if (replication == 1L || keep_tests) tests else NULL)
  }, error = function(e) {
    stop(sprintf("EEG replication %d (seed %d) failed: %s. No partial average was returned.",
                 replication, seeds[replication], conditionMessage(e)), call. = FALSE)
  })
}

# R is the replication count. Seeds are seed, seed+1, ..., seed+R-1.
# Changing R or n_cores preserves the data generated for each replication.
# n_cores=1 is sequential; n_cores>1 uses the parallel package included with R.
# No external files are sourced and no files are written by this function.
run_EEG_simulated = function(R = 1L, m = 50000L, burnin = 1000L,
                             TV_size = 3, sigma0 = 25, mean_level = 0,
                             seed = 20261015L, alpha = 0.05, nu = 2L, L = 1L,
                             Kernel = "Kernel_Triangular", n_cores = 1L,
                             keep_tests = FALSE, verbose = TRUE,
                             onset = 0.68, rise_width = 0.08,
                             fall_start = 0.88, fall_width = 0.06,
                             late_fraction = 0.5) {
  .eeg_check_integer(R, "R", 1L)
  .eeg_check_integer(m, "m", 20L)
  .eeg_check_integer(burnin, "burnin", 0L)
  .eeg_check_integer(seed, "seed", 0L)
  .eeg_check_integer(n_cores, "n_cores", 1L)
  .eeg_check_integer(nu, "nu", 1L)
  .eeg_check_integer(L, "L", 0L)
  .eeg_check_flag(keep_tests, "keep_tests")
  .eeg_check_flag(verbose, "verbose")
  if (seed + R - 1 > .Machine$integer.max) stop("seed + R - 1 exceeds R's seed limit.")
  if (m + burnin > .Machine$integer.max) stop("m + burnin is too large.")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must be strictly between zero and one.")
  if (m %% 2 != 0 || m < 16 * (12 * (2 * nu + 1) - 1) + 4)
    stop("m must be even and large enough for the selected nu (at least 948 for nu=2).")
  seeds = as.integer(seed + seq_len(R) - 1)
  simulation_args = list(m = m, burnin = burnin, TV_size = TV_size,
                         sigma0 = sigma0, mean_level = mean_level,
                         onset = onset, rise_width = rise_width,
                         fall_start = fall_start, fall_width = fall_width,
                         late_fraction = late_fraction, return_details = TRUE)
  fit_args = list(nu = nu, L = L, Kernel = Kernel)
  workers = min(as.integer(n_cores), as.integer(R))
  start = proc.time()[["elapsed"]]
  if (workers == 1L) {
    runs = lapply(seq_len(R), .eeg_run_one, simulation_args = simulation_args,
                  fit_args = fit_args, seeds = seeds, alpha = alpha,
                  keep_tests = keep_tests, verbose = verbose)
  } else {
    if (verbose) message(sprintf("Running %d replications on %d workers; waiting for all results...", R, workers))
    cl = parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    runs = parallel::parLapply(cl, seq_len(R), .eeg_run_one,
                              simulation_args = simulation_args, fit_args = fit_args,
                              seeds = seeds, alpha = alpha, keep_tests = keep_tests,
                              verbose = FALSE)
  }
  summary = summarize_EEG_p12(lapply(runs, `[[`, "score"))
  summary$per_replication_metrics$seed = seeds
  summary$per_replication_metrics$elapsed_seconds = vapply(runs, `[[`, numeric(1), "elapsed_seconds")
  summary$replication_diagnostics = data.frame(
    replication = seq_len(R), seed = seeds,
    frequencies = vapply(runs, function(z) paste(z$frequencies, collapse = ","), character(1)),
    bandwidths = vapply(runs, function(z) paste(z$bandwidths, collapse = ","), character(1)))
  summary$scores = lapply(runs, `[[`, "score")
  summary$first_simulation = runs[[1]]$simulation
  summary$first_tests = runs[[1]]$tests
  summary$all_tests = if (keep_tests) lapply(runs, `[[`, "tests") else NULL
  summary$truth = runs[[1]]$simulation$truth
  summary$settings = c(list(R = R, seed = seed, seeds = seeds, alpha = alpha,
                           n_cores = workers, keep_tests = keep_tests),
                      simulation_args, fit_args)
  summary$elapsed_seconds = proc.time()[["elapsed"]] - start
  if (verbose) message(sprintf("Finished all %d replications in %.1f seconds.", R, summary$elapsed_seconds))
  summary
}

list(simulate = sim.tvVAR_p12_eeg, spec = p12_eeg_spec,
     precision = p12_precision_blocks, fit = fit_EEG_simulated,
     score = score_EEG_p12, summarize = summarize_EEG_p12,
     run = run_EEG_simulated)
}, envir = new.env(parent = asNamespace("stats")))

# Public functions; the copied legacy estimator remains inside its private environment.
sim.tvVAR_p12_eeg = .EEG_simulated_module$simulate
p12_eeg_spec = .EEG_simulated_module$spec
p12_precision_blocks = .EEG_simulated_module$precision
fit_EEG_simulated = .EEG_simulated_module$fit
score_EEG_p12 = .EEG_simulated_module$score
summarize_EEG_p12 = .EEG_simulated_module$summarize
run_EEG_simulated = .EEG_simulated_module$run

# -----------------------------------------------------------------------------
# RUN THE STUDY: edit R below to choose the number of replications.
# Default seed reproduces the previously generated p12_seed_20261015 example.
# For functions only: options(NonStGM.EEG_simulated.run_example = FALSE) before source().
# -----------------------------------------------------------------------------
if (isTRUE(getOption("NonStGM.EEG_simulated.run_example", interactive()))) {
  EEG_simulated_results = run_EEG_simulated(
    R = 1L,                 # Set to 100L, for example, to average 100 replications.
    m = 50000L,
    burnin = 1000L,
    TV_size = 3,
    sigma0 = 25,
    seed = 20261015L,
    alpha = 0.05,
    nu = 2L,
    L = 1L,
    n_cores = 1L,           # Optionally use 2L or 4L when R > 1.
    keep_tests = FALSE,     # TRUE retains every replication's raw test table.
    verbose = TRUE
  )

  # These always refer to the FIRST replication, even when R > 1.
  EEG_simulated_x = EEG_simulated_results$first_simulation$x
  EEG_simulated_tests = EEG_simulated_results$first_tests
  print(EEG_simulated_results$truth$edges)
  print(EEG_simulated_results$truth$nodes)
  print(EEG_simulated_results$metric_summary, row.names = FALSE)

  if (interactive()) {
    graphics::plot.ts(EEG_simulated_x[, 1], xlab = "Time", ylab = "X1",
                      main = "Simulated EEG-like VAR: first replication")
    graphics::abline(v = 34000, col = "red", lty = 2)
  }
}

# Useful results after the example:
# EEG_simulated_results$per_replication_metrics  # one row per replication
# EEG_simulated_results$metric_summary           # mean, MCSE, n_valid per metric
# EEG_simulated_results$edge_detection_frequency # all 66 pairs, detection percentages
# EEG_simulated_results$node_detection_frequency # all 12 node labels
# EEG_simulated_results$scores[[1]]               # decisions for replication 1
# EEG_simulated_results$replication_diagnostics  # chosen frequencies and bandwidths
# EEG_simulated_results$truth                     # actual graph and node truth
#
# To run more after sourcing without repeating the example:
# study_100 = run_EEG_simulated(R = 100L, n_cores = 4L)
# study_100$metric_summary
#
# To fit an already-generated example with the same one-step method:
# tests = fit_EEG_simulated(EEG_simulated_x)
# score = score_EEG_p12(tests, EEG_simulated_results$truth, nu = 2L,
#                       frequencies = attr(tests, "selected_frequencies"))
