# ============================================================
# p=5 VAR SIMULATION FUNCTION
# ============================================================
p =5
sim.tvVAR_p5 <- function(burnin, m, TV_size) {

  A <- matrix(c(
    0.5,  0,    0,    0,    0,
    0.2,  0.6,  0,    0,    0,
    0,    0.3,  0.5,  0,    0,
    0,    0,    0.25, 0.55, 0,
    0.15, 0,    0,    0.2,  0.5
  ), ncol = 5, byrow = TRUE)

  n <- m + burnin
  p <- 5
  x <- matrix(rnorm(n * p), ncol = p)
  x1 <- x

  # Time-varying coefficient for A[1,1] - same formula as p=3
  st <- 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)

  for (tt in 2:n) {
    A.t <- A
    A.t[1, 1] <- st[tt]  # Only node 1 is time-varying
    temp <- A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] <- c(temp)
  }

  x2 <- x1[-c(1:burnin), ]
  return(x2)
}

# ============================================================
# GRAPH RECOVERY FUNCTION FOR p=5
# ============================================================
# True edges:
#   Time-varying: (1,1)
#   Off-diagonal: (1,2), (2,3), (3,4), (4,5), (1,5)
#   Total: 6 true edges, 5 off-diagonal

graph_recovery_p5 <- function(Test_tibble, alpha = 0.05) {

  # Filter to lower triangular
  Test_tibble <- Test_tibble |>
    filter(a <= c)

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
  # Step 1: Within each (a, c, r, k, i) group
  Test_tibble_twostep <- Test_tibble |>
    group_by(a, c, r, k, i) |>
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
      global_detect_11_tv = any(a == 1 & c == 1 & r > 0 & global_BY_min < alpha),
      global_detect_12 = any(a == 1 & c == 2 & global_BY_min < alpha) |
        any(a == 2 & c == 1 & global_BY_min < alpha),
      global_detect_23 = any(a == 2 & c == 3 & global_BY_min < alpha) |
        any(a == 3 & c == 2 & global_BY_min < alpha),
      global_detect_34 = any(a == 3 & c == 4 & global_BY_min < alpha) |
        any(a == 4 & c == 3 & global_BY_min < alpha),
      global_detect_45 = any(a == 4 & c == 5 & global_BY_min < alpha) |
        any(a == 5 & c == 4 & global_BY_min < alpha),
      global_detect_15 = any(a == 1 & c == 5 & global_BY_min < alpha) |
        any(a == 5 & c == 1 & global_BY_min < alpha),

      # Count number of significant off-diagonal edges
      global_n_edges = nrow(distinct(data.frame(
        a2 = pmin(a[a != c & global_BY_min < alpha],
                  c[a != c & global_BY_min < alpha]),
        c2 = pmax(a[a != c & global_BY_min < alpha],
                  c[a != c & global_BY_min < alpha])
      ))),

      # For Two-Step BY method - detect each true edge
      twostep_detect_11_tv = any(a == 1 & c == 1 & r > 0 & twostep_BY < alpha),
      twostep_detect_12 = any(a == 1 & c == 2 & twostep_BY < alpha) |
        any(a == 2 & c == 1 & twostep_BY < alpha),
      twostep_detect_23 = any(a == 2 & c == 3 & twostep_BY < alpha) |
        any(a == 3 & c == 2 & twostep_BY < alpha),
      twostep_detect_34 = any(a == 3 & c == 4 & twostep_BY < alpha) |
        any(a == 4 & c == 3 & twostep_BY < alpha),
      twostep_detect_45 = any(a == 4 & c == 5 & twostep_BY < alpha) |
        any(a == 5 & c == 4 & twostep_BY < alpha),
      twostep_detect_15 = any(a == 1 & c == 5 & twostep_BY < alpha) |
        any(a == 5 & c == 1 & twostep_BY < alpha),

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
      # Global BY: Perfect recovery (all 6 true edges, exactly 5 off-diagonal)
      global_all_true_detected = global_detect_11_tv & global_detect_12 &
        global_detect_23 & global_detect_34 &
        global_detect_45 & global_detect_15,
      global_perfect = global_all_true_detected & (global_n_edges == 5),

      # Two-Step BY: Perfect recovery
      twostep_all_true_detected = twostep_detect_11_tv & twostep_detect_12 &
        twostep_detect_23 & twostep_detect_34 &
        twostep_detect_45 & twostep_detect_15,
      twostep_perfect = twostep_all_true_detected & (twostep_n_edges == 5),

      # Sensitivity: proportion of 6 true edges detected
      global_sensitivity = (global_detect_11_tv + global_detect_12 + global_detect_23 +
                              global_detect_34 + global_detect_45 + global_detect_15) / 6,
      twostep_sensitivity = (twostep_detect_11_tv + twostep_detect_12 + twostep_detect_23 +
                               twostep_detect_34 + twostep_detect_45 + twostep_detect_15) / 6
    )

  # Overall accuracy metrics
  overall_accuracy <- tibble(
    Method = c("Global BY", "Two-Step BY"),

    # Perfect graph recovery rate
    Perfect_Recovery = c(
      mean(accuracy_summary$global_perfect),
      mean(accuracy_summary$twostep_perfect)
    ),

    # All true edges detected (may have false positives)
    All_True_Detected = c(
      mean(accuracy_summary$global_all_true_detected),
      mean(accuracy_summary$twostep_all_true_detected)
    ),

    # Individual edge detection rates
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
    Detect_15 = c(
      mean(accuracy_summary$global_detect_15),
      mean(accuracy_summary$twostep_detect_15)
    ),

    # Sensitivity (proportion of true edges detected)
    Mean_Sensitivity = c(
      mean(accuracy_summary$global_sensitivity),
      mean(accuracy_summary$twostep_sensitivity)
    ),

    # Mean number of detected off-diagonal edges
    Mean_N_Edges = c(
      mean(accuracy_summary$global_n_edges),
      mean(accuracy_summary$twostep_n_edges)
    ),

    # False positives (extra off-diagonal edges beyond 5 true ones)
    Mean_False_Positives = c(
      mean(pmax(0, accuracy_summary$global_n_edges - 5)),
      mean(pmax(0, accuracy_summary$twostep_n_edges - 5))
    )
  )

  return(overall_accuracy)
}

# ============================================================
# SIMULATION RUNNER FOR p=5
# ============================================================

run_NonStGM_simulation_p5 <- function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                      TV_size = 0.6, nu = 2, L = 1,
                                      Kernel = 'Kernel_Triangular', seed = 0) {

  # Collect all Test_tibbles across replications
  All_Test_tibble <- NULL

  set.seed(seed)

  for (iter in 1:R) {

    # Generate p=5 data
    x <- sim.tvVAR_p5(burnin, m, TV_size)

    # Get Test_tibble from NonStGM_adj
    tib <- NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)

    # Update iteration index to reflect replication number
    tib <- tib |>
      mutate(i = iter)

    # Combine with all replications
    All_Test_tibble <- rbind(All_Test_tibble, tib)

    if (iter %% 20 == 0) cat("Completed iteration:", iter, "/", R, "\n")
  }

  # Apply graph_recovery_p5 to combined tibble
  results <- graph_recovery_p5(All_Test_tibble, alpha = alpha)

  return(results)
}

# ============================================================
# RUN THE SIMULATION
# ============================================================

results_p5 <- run_NonStGM_simulation_p5(
  R = 20,
  alpha = 0.5,
  burnin = 20,
  m = 4096,
  TV_size = 0.6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 0
)

print(results_p5)
