# ============================================================
# p=7 VAR SIMULATION FUNCTION
# ============================================================
source('extractVariables.R')
source('exampleVAR.R')
source('KernelWeights.R')
source('Combined_MK_Estimation.R')
source('DFTransform.R')
source('R_Hat_Creation.R')
source('BetaFunctions.R')
source('VarianceFunctions.R')
#source('NonStGM_adj.R')
p <- 7

sim.tvVAR_p7 <- function(burnin, m, TV_size) {
  # Circular structure: 1-2-3-4-5-6-7-1
  # All diagonal = 0.5, all off-diagonal coefficients = 0.3
  A <- matrix(c(
    0.5,  0,    0,    0,    0,    0,    0,
    0.3,  0.5,  0,    0,    0,    0,    0,
    0,    0.3,  0.5,  0,    0,    0,    0,
    0,    0,    0.3,  0.5,  0,    0,    0,
    0,    0,    0,    0.3,  0.5,  0,    0,
    0,    0,    0,    0,    0.3,  0.5,  0,
    0.3,  0,    0,    0,    0,    0.3,  0.5
  ), ncol = 7, byrow = TRUE)

  n <- m + burnin
  p <- 7
  x <- matrix(rnorm(n * p), ncol = p)
  x1 <- x

  # Time-varying coefficient for A[1][1] - same formula as p=3 and p=5
  st <- 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)

  for (tt in 2:n) {
    A.t <- A
    A.t[1][1] <- st[tt]  # Only node 1 is time-varying
    temp <- A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] <- c(temp)
  }

  x2 <- x1[-c(1:burnin), ]
  return(x2)
}

check_tvVAR_stability_p7 <- function(burnin, m, TV_size) {
  A <- matrix(c(
    0.5,  0,    0,    0,    0,    0,    0,
    0.3,  0.5,  0,    0,    0,    0,    0,
    0,    0.3,  0.5,  0,    0,    0,    0,
    0,    0,    0.3,  0.5,  0,    0,    0,
    0,    0,    0,    0.3,  0.5,  0,    0,
    0,    0,    0,    0,    0.3,  0.5,  0,
    0.3,  0,    0,    0,    0,    0.3,  0.5
  ), ncol = 7, byrow = TRUE)

  n <- burnin + m
  st <- 0.3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  max_mod <- sapply(st, function(s) { A[1][1] <- s; max(Mod(eigen(A)$values)) })

  cat("A[1][1] range:", round(range(st), 3), "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}

# Usage
# check_tvVAR_stability_p7(20, 4096, 0.6)

# ============================================================
# GRAPH RECOVERY FUNCTION FOR p=7
# ============================================================
# True edges:
#   Time-varying: (1,1)
#   Off-diagonal: (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (1,7)
#   Total: 8 true edges, 7 off-diagonal

graph_recovery_p7 <- function(Test_tibble, alpha = 0.05) {
  # Filter to lower triangular
  Test_tibble <- Test_tibble |>
    filter(a <= c)

  # Filter r values: r = 0 for off-diagonal (a != c), r = 1 for diagonal (a == c)
  Test_tibble <- Test_tibble |>
    filter((a != c & r == 0) | (a == c & r == 1))

  # Global BY adjustment
  Test_tibble_global <- Test_tibble |>
    group_by(i) |>
    group_modify(\(df, key) {
      all_pvals <- c(df$Re, df$Im)
      all_adj <- p.adjust(all_pvals, method = "BY")
      n_tests <- nrow(df)
      df$Re_global_BY <- all_adj[1:n_tests]
      df$Im_global_BY <- all_adj[(n_tests + 1):(2 * n_tests)]
      df$global_BY_min <- pmin(df$Re_global_BY, df$Im_global_BY)
      df
    }) |>
    ungroup()

  # Two-step BY adjustment
  # Step 1: Within each (a, c, k, i) group
  Test_tibble_twostep <- Test_tibble |>
    group_by(a, c, k, i) |>
    group_modify(\(df, key) {
      all_pvals <- c(df$Re, df$Im)
      all_adj <- p.adjust(all_pvals, method = "BY")
      n_tests <- nrow(df)
      df$Re_step1_BY <- all_adj[1:n_tests]
      df$Im_step1_BY <- all_adj[(n_tests + 1):(2 * n_tests)]
      df$step1_min <- pmin(df$Re_step1_BY, df$Im_step1_BY)
      df
    }) |>
    ungroup()

  # Step 2: Apply BY correction across all step1_min values
  Test_tibble_twostep <- Test_tibble_twostep |>
    mutate(twostep_BY = p.adjust(step1_min, method = "BY"))

  # Step 3: Merge results
  Test_tibble_final <- Test_tibble |>
    left_join(
      Test_tibble_global |> select(a, c, r, k, i, global_BY_min),
      by = c("a", "c", "r", "k", "i")
    ) |>
    left_join(
      Test_tibble_twostep |> select(a, c, r, k, i, twostep_BY),
      by = c("a", "c", "r", "k", "i")
    )

  # Graph recovery metrics
  graph_recovery_results <- Test_tibble_final |>
    group_by(i) |>
    summarise(
      # For Global BY method - detect each true edge
      global_detect_11_tv = any(a == 1 & c == 1 & global_BY_min < alpha),
      global_detect_12 = any(a == 1 & c == 2 & global_BY_min < alpha) |
        any(a == 2 & c == 1 & global_BY_min < alpha),
      global_detect_23 = any(a == 2 & c == 3 & global_BY_min < alpha) |
        any(a == 3 & c == 2 & global_BY_min < alpha),
      global_detect_34 = any(a == 3 & c == 4 & global_BY_min < alpha) |
        any(a == 4 & c == 3 & global_BY_min < alpha),
      global_detect_45 = any(a == 4 & c == 5 & global_BY_min < alpha) |
        any(a == 5 & c == 4 & global_BY_min < alpha),
      global_detect_56 = any(a == 5 & c == 6 & global_BY_min < alpha) |
        any(a == 6 & c == 5 & global_BY_min < alpha),
      global_detect_67 = any(a == 6 & c == 7 & global_BY_min < alpha) |
        any(a == 7 & c == 6 & global_BY_min < alpha),
      global_detect_17 = any(a == 1 & c == 7 & global_BY_min < alpha) |
        any(a == 7 & c == 1 & global_BY_min < alpha),
      # Count number of significant off-diagonal edges
      global_n_edges = nrow(distinct(data.frame(
        a2 = pmin(a[a != c & global_BY_min < alpha],
                  c[a != c & global_BY_min < alpha]),
        c2 = pmax(a[a != c & global_BY_min < alpha],
                  c[a != c & global_BY_min < alpha])
      ))),
      # For Two-Step BY method - detect each true edge
      twostep_detect_11_tv = any(a == 1 & c == 1 & twostep_BY < alpha),
      twostep_detect_12 = any(a == 1 & c == 2 & twostep_BY < alpha) |
        any(a == 2 & c == 1 & twostep_BY < alpha),
      twostep_detect_23 = any(a == 2 & c == 3 & twostep_BY < alpha) |
        any(a == 3 & c == 2 & twostep_BY < alpha),
      twostep_detect_34 = any(a == 3 & c == 4 & twostep_BY < alpha) |
        any(a == 4 & c == 3 & twostep_BY < alpha),
      twostep_detect_45 = any(a == 4 & c == 5 & twostep_BY < alpha) |
        any(a == 5 & c == 4 & twostep_BY < alpha),
      twostep_detect_56 = any(a == 5 & c == 6 & twostep_BY < alpha) |
        any(a == 6 & c == 5 & twostep_BY < alpha),
      twostep_detect_67 = any(a == 6 & c == 7 & twostep_BY < alpha) |
        any(a == 7 & c == 6 & twostep_BY < alpha),
      twostep_detect_17 = any(a == 1 & c == 7 & twostep_BY < alpha) |
        any(a == 7 & c == 1 & twostep_BY < alpha),
      # Count number of significant off-diagonal edges for two-step
      twostep_n_edges = nrow(distinct(data.frame(
        a2 = pmin(a[a != c & twostep_BY < alpha],
                  c[a != c & twostep_BY < alpha]),
        c2 = pmax(a[a != c & twostep_BY < alpha],
                  c[a != c & twostep_BY < alpha])
      ))),
      .groups = "drop"
    )

  # Accuracy summary
  accuracy_summary <- graph_recovery_results |>
    mutate(
      # Global BY: Perfect recovery (all 8 true edges, exactly 7 off-diagonal) ( removed global_detect_11_tv )
      global_all_true_detected =  global_detect_12 &
        global_detect_23 & global_detect_34 &
        global_detect_45 & global_detect_56 &
        global_detect_67 & global_detect_17,
      global_perfect = global_all_true_detected & (global_n_edges == 7),
      # Two-Step BY: Perfect recovery
      twostep_all_true_detected = twostep_detect_11_tv & twostep_detect_12 &
        twostep_detect_23 & twostep_detect_34 &
        twostep_detect_45 & twostep_detect_56 &
        twostep_detect_67 & twostep_detect_17,
      twostep_perfect = twostep_all_true_detected & (twostep_n_edges == 7),
      # Sensitivity: proportion of 8 true edges detected
      global_sensitivity = (global_detect_11_tv + global_detect_12 + global_detect_23 +
                              global_detect_34 + global_detect_45 + global_detect_56 +
                              global_detect_67 + global_detect_17) / 8,
      twostep_sensitivity = (twostep_detect_11_tv + twostep_detect_12 + twostep_detect_23 +
                               twostep_detect_34 + twostep_detect_45 + twostep_detect_56 +
                               twostep_detect_67 + twostep_detect_17) / 8
    )

  # Overall accuracy metrics
  overall_accuracy <- tibble(
    Method = c("Global BY", "Two-Step BY"),
    Perfect_Recovery = c(
      mean(accuracy_summary$global_perfect),
      mean(accuracy_summary$twostep_perfect)
    ),
    All_True_Detected = c(
      mean(accuracy_summary$global_all_true_detected),
      mean(accuracy_summary$twostep_all_true_detected)
    ),
    Detect_11_TV = c(
      mean(accuracy_summary$global_detect_11_tv),
      mean(accuracy_summary$twostep_detect_11_tv)
    ),
    Detect_12 = c(
      mean(accuracy_summary$global_detect_12),
      mean(accuracy_summary$twostep_detect_12)
    ),
    Detect_23 = c(
      mean(accuracy_summary$global_detect_23),
      mean(accuracy_summary$twostep_detect_23)
    ),
    Detect_34 = c(
      mean(accuracy_summary$global_detect_34),
      mean(accuracy_summary$twostep_detect_34)
    ),
    Detect_45 = c(
      mean(accuracy_summary$global_detect_45),
      mean(accuracy_summary$twostep_detect_45)
    ),
    Detect_56 = c(
      mean(accuracy_summary$global_detect_56),
      mean(accuracy_summary$twostep_detect_56)
    ),
    Detect_67 = c(
      mean(accuracy_summary$global_detect_67),
      mean(accuracy_summary$twostep_detect_67)
    ),
    Detect_17 = c(
      mean(accuracy_summary$global_detect_17),
      mean(accuracy_summary$twostep_detect_17)
    ),
    Mean_Sensitivity = c(
      mean(accuracy_summary$global_sensitivity),
      mean(accuracy_summary$twostep_sensitivity)
    ),
    Mean_N_Edges = c(
      mean(accuracy_summary$global_n_edges),
      mean(accuracy_summary$twostep_n_edges)
    ),
    Mean_False_Positives = c(
      mean(pmax(0, accuracy_summary$global_n_edges - 7)),
      mean(pmax(0, accuracy_summary$twostep_n_edges - 7))
    )
  )

  return(overall_accuracy)
}

# ============================================================
# SIMULATION RUNNER FOR p=7
# ============================================================
run_NonStGM_simulation_p7 <- function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                      TV_size = 0.6, nu = 2, L = 1,
                                      Kernel = 'Kernel_Triangular', seed = 0) {
  # Collect all Test_tibbles across replications
  All_Test_tibble <- NULL
  set.seed(seed)

  for (iter in 1:R) {
    # Generate p=7 data
    x <- sim.tvVAR_p7(burnin, m, TV_size)

    # Get Test_tibble from NonStGM_adj
    tib <- NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)

    # Update iteration index to reflect replication number
    tib <- tib |>
      mutate(i = iter)

    # Combine with all replications
    All_Test_tibble <- rbind(All_Test_tibble, tib)
  }
  # Apply graph_recovery_p7 to combined tibble
  results <- graph_recovery_p7(All_Test_tibble, alpha = alpha)

  return(list(results,All_Test_tibble))
}

# ============================================================
# RUN THE SIMULATION
# ============================================================
results_p7 <- run_NonStGM_simulation_p7(
  R = 100,
  alpha = 0.05,
  burnin = 20,
  m = 2^12,
  TV_size = 0.6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 0
)

print(results_p7[[1]])
hi = results_p7[[2]]
