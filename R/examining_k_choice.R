# Define K range bins to examine
k_ranges <- list(
  "all"      = c(-Inf, Inf),
  "low"      = c(1400, 1500),
  "mid"      = c(1501, 1600),
  "high"     = c(1601, 1700),
  "very_high"= c(1700, max(results_p7[[2]]$k, na.rm = TRUE))
)

filtered_results <- results_p7[[2]] |>
  filter(a == 1 & c == 1 & r == 1)

# Compute p_min for each k range
p_value_by_k_range <- bind_rows(lapply(names(k_ranges), function(rng) {
  k_lo <- k_ranges[[rng]][1]
  k_hi <- k_ranges[[rng]][2]

  filtered_results |>
    filter(k >= k_lo & k <= k_hi) |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      n_k   = n(),
      .groups = "drop"
    ) |>
    mutate(k_range = rng)
}))

# Faceted histogram
library(ggplot2)

ggplot(p_value_by_k_range, aes(x = p_min)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  facet_wrap(~ k_range, scales = "free_y") +
  labs(
    title = "P-value Distribution by K Range (a=1, c=1, r=1)",
    subtitle = "Checking if signal is sparse across frequencies",
    x = "BY-adjusted p-value (min over Re & Im)",
    y = "Count"
  ) +
  theme_minimal()

# Summary table: detection rate per k range
detection_summary <- p_value_by_k_range |>
  group_by(k_range) |>
  summarise(
    n_nodes      = n(),
    n_significant = sum(p_min < 0.05),
    detection_rate = mean(p_min < 0.05),
    median_p     = median(p_min),
    mean_p       = mean(p_min),
    .groups = "drop"
  ) |>
  arrange(detection_rate)

print(detection_summary)

per_k_signal <- filtered_results |>
  mutate(k_bin = floor(k / 10) * 10) |>  # bucket into groups of 10
  group_by(k_bin) |>
  summarise(
    mean_neg_log_p = mean(-log10(pmin(Re, Im)), na.rm = TRUE),
    median_neg_log_p = median(-log10(pmin(Re, Im)), na.rm = TRUE),
    min_p = min(pmin(Re, Im), na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

ggplot(per_k_signal, aes(x = k_bin, y = mean_neg_log_p)) +
  geom_col(fill = "steelblue", width = 8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = max(per_k_signal$k_bin) * 0.95, y = -log10(0.05) + 0.15,
           label = "p = 0.05", color = "red", size = 3.5) +
  labs(
    title = "Average Signal Strength by Frequency Bin (a=1, c=1, r=1)",
    subtitle = "Frequencies bucketed in groups of 10, averaged over all nodes",
    x = "Frequency bin (k)",
    y = "Mean -log10(p)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Histogram of realized k frequencies, bucketed by 10
k_hist <- filtered_results |>
  mutate(k_bin = floor(k / 10) * 10) |>
  group_by(k_bin) |>
  summarise(count = n(), .groups = "drop")

# Plot: histogram of selected k values
ggplot(k_hist, aes(x = k_bin, y = count)) +
  geom_col(fill = "darkorange", width = 8, alpha = 0.8) +
  labs(
    title = "Distribution of Realized K Frequencies (a=1, c=1, r=1)",
    subtitle = "How often each frequency bin was selected by K selection procedure",
    x = "Frequency bin (k)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combined dual-panel plot for direct comparison
library(patchwork)

p1 <- ggplot(k_hist, aes(x = k_bin, y = count)) +
  geom_col(fill = "darkorange", width = 8, alpha = 0.8) +
  labs(title = "K Selection Frequency", x = "Frequency bin (k)", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(per_k_signal, aes(x = k_bin, y = mean_neg_log_p)) +
  geom_col(fill = "steelblue", width = 8, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Average Signal Strength", x = "Frequency bin (k)", y = "Mean -log10(p)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 / p2 + plot_annotation(
  title = "K Selection vs Signal Strength (a=1, c=1, r=1)",
  subtitle = "Frequencies bucketed in groups of 10"
)
