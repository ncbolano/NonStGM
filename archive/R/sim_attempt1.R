
# True edges for NonStGM (detects NON-STATIONARITY):
# Only (1,1) is time-varying, connections (1,2) and (2,3) exist
true_edges = tibble(
  node1 = c(1, 1, 2),
  node2 = c(1, 2, 3)
)
n_true_edges = nrow(true_edges)

# Total possible edges (including self-loops) for p=3:
# (1,1), (1,2), (1,3), (2,2), (2,3), (3,3) = 6 edges
p = 3
total_possible_edges = p * (p + 1) / 2  # 6 for p=3

# Simulation parameters
n_sim = 100
alpha = 0.05
burnin = 500
m = 4096
TV_size = 0.6

# Pre-allocate results
fdr_vec = numeric(n_sim)
power_vec = numeric(n_sim)
perfect_vec = logical(n_sim)
n_detected_vec = numeric(n_sim)
tp_vec = numeric(n_sim)
fp_vec = numeric(n_sim)

set.seed(0)

for (iter in 1:n_sim) {

  x = sim.tvVAR(burnin, m, TV_size)
  significant_edges = NonStGM(x, "Kernel_Triangular", 2, 1, alpha = alpha)

  n_detected = nrow(significant_edges)
  n_detected_vec[iter] = n_detected

  if (n_detected == 0) {
    tp = 0
    fp = 0
    fdr_vec[iter] = 0
    power_vec[iter] = 0
    perfect_vec[iter] = FALSE
  } else {
    # Create edge keys for comparison
    detected_keys = paste(significant_edges$node1, significant_edges$node2, sep = "_")
    true_keys = paste(true_edges$node1, true_edges$node2, sep = "_")

    tp = sum(detected_keys %in% true_keys)

    fp = n_detected - tp

    tp_vec[iter] = tp
    fp_vec[iter] = fp

    fdr_vec[iter] = fp / n_detected

    power_vec[iter] = tp / n_true_edges

    perfect_vec[iter] = setequal(detected_keys, true_keys)
  }

  if (iter %% 100 == 0) cat("Completed iteration:", iter, "\n")
}

# Summary results
summary_results <- tibble(
  Metric = c("Number of Variables (p)",
             "Number of Observations (n)",
             "Significance Level (alpha)",
             "Number of Simulations",
             "Number of True Edges",
             "Avg Detected Edges",
             "Avg True Positives",
             "Avg False Positives",
             "Average FDR",
             "Average Power",
             "Perfect Recovery Rate"),
  Value = c(p,
            m,
            alpha,
            n_sim,
            n_true_edges,
            round(mean(n_detected_vec), 3),
            round(mean(tp_vec), 3),
            round(mean(fp_vec), 3),
            round(mean(fdr_vec), 4),
            round(mean(power_vec), 4),
            round(mean(perfect_vec), 4))
)

print(summary_results)
