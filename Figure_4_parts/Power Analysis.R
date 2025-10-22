# R Code for Power Analysis in Review Paper
# Sarah Ruckman

# Load required library
library(ggplot2)
library(dplyr)

#  ===== The simulation function  =====
calculate_power <- function(sample_size, pve, num_tests, num_simulations = 1000) {
  
  # Step 1: Calculate the Bonferroni-corrected significance threshold
  alpha <- 0.05 / num_tests
  
  # Step 2: Define model parameters
  total_variance <- 100
  effect_size_a <- sqrt(2 * pve * total_variance)
  environmental_sd <- sqrt((1 - pve) * total_variance)
  
  # Step 3: Initialize counter for significant results
  significant_results_count <- 0
  
  # Step 4: Run simulation loop
  for (i in 1:num_simulations) {
    
    # a. Simulate Genotypes (G)
    G <- rbinom(n = sample_size, size = 2, prob = 0.5)
    
    # b. Simulate Phenotypes (Y)
    E <- rnorm(n = sample_size, mean = 0, sd = environmental_sd)
    Y <- (G - 1) * effect_size_a + E
    
    # c. Perform Statistical Test
    model <- lm(Y ~ G)
    p_value <- summary(model)$coefficients[2, 4]
    
    # d. Check for significance
    if (p_value < alpha) {
      significant_results_count <- significant_results_count + 1
    }
  }
  
  # Step 5: Calculate power
  power <- significant_results_count / num_simulations
  
  # Step 6: Return result
  return(power)
}

# ===== MAIN ANALYSIS =====

# Define the scenarios from your table
scenarios <- data.frame(
  Line = c("MPP", "MPP", "MPP", "MPP", "Outbred", "Outbred", "Outbred", "Outbred"),
  Phenotypes = c(1, 1, 50000, 50000, 1, 1, 50000, 50000),
  Analysis = c("WG", "Cis-only", "WG", "Cis-only", "WG", "Cis-only", "WG", "Cis-only"),
  num_tests = c(6000, 2, 3e8, 1e5, 2.4e6, 800, 1.2e11, 4e7),
  stringsAsFactors = FALSE
)

# Create a readable label for each scenario
scenarios$scenario_name <- paste(scenarios$Line, scenarios$Phenotypes, scenarios$Analysis, sep=" ")

# Define the parameter ranges
pve_values <- c(0.01, 0.02, 0.05, 0.10)  # 2%, 5%, 10% variance explained
sample_sizes <- c(100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 5000)

# Create an empty data frame to store all results
results_df <- data.frame()

# Loop over each scenario (each row in your table)
for (i in 1:nrow(scenarios)) {
  
  # Loop over each Percent Variance Explained (PVE)
  for (pve in pve_values) {
    
    # Loop over each Sample Size
    for (n in sample_sizes) {
      
      # Print progress (helpful for long simulations)
      cat(sprintf("Running: %s | PVE: %.2f | N: %d\n", 
                  scenarios$scenario_name[i], pve, n))
      
      # Run the simulation
      current_power <- calculate_power(
        sample_size = n,
        pve = pve,
        num_tests = scenarios$num_tests[i],
        num_simulations = 1000  # You can adjust this
      )
      
      # Store the result
      temp_result <- data.frame(
        pve = pve,
        sample_size = n,
        scenario = scenarios$scenario_name[i],
        power = current_power,
        stringsAsFactors = FALSE
      )
      
      results_df <- rbind(results_df, temp_result)
    }
  }
}

# Save the results to a CSV file
write.csv(results_df, "power_analysis_results.csv", row.names = FALSE)

# ===== CREATE THE PLOTS =====
results_df <- read.csv("C:\\Users\\sarah\\OneDrive - UC Irvine\\Long Lab\\Short Read Review\\power_analysis_5000.csv", header = T)


# Create labels for the facets
pve_labels <- c(
  "0.02" = "2% Variance Explained",
  "0.05" = "5% Variance Explained",
  "0.1" = "10% Variance Explained"
)

# CHOOSE YOUR 8 COLORS HERE:
my_colors <- c(
  "#00008B",  # Color 1 = MPP 1 WG
  "#CD2626",  # Color 3 = MPP 50k Cis
  "#0000EE",  # Color 2 = Outbred 1 WG
  "palevioletred1",  # Color 5 = Outbred 50k Cis
  "#1C86EE",  # Color 4 = MPP 50k WG
  "#87CEFA",  # Color 6 = Outbred 50k WG
  "#6E8B3D",  # Color 7 
  "#FFC0CB"   # Color 8
)

my_colors <- c(
  "#1C86EE",  # Color 2 - Dodger Blue 2 = MPP 1 WG
  "#BF3EFF",  # Color 1 - Dark orchid 1   = MPP 50k Cis
  "#1C86EE",  # Color 2 - Dodger Blue 2  = Outbred 1 WG
  "#BF3EFF",  # Color 1 - Dark orchid 1 = Outbred 50k Cis
  "#1C86EE",  # Color 2 - Dodger Blue 2 = MPP 50k WG
  "#1C86EE",  # Color 2 - Dodger Blue 2 = Outbred 50k WG
  "#BF3EFF",  # Color 1 - Dark orchid 1
  "#1C86EE"   # Color 2 - Dodger Blue 2 
)
# FILTER OUT the non-sensical cis-only single phenotype designs AND 1% variance
results_df_filtered <- subset(results_df, 
                              !(scenario %in% c("MPP 1 Cis-only", "Outbred 1 Cis-only")) & pve != 0.01)

# Add a column to identify the Line type for shapes
results_df_filtered$line_type <- ifelse(grepl("^MPP", results_df_filtered$scenario), "MPP", "Outbred")
results_df_filtered$phenotype_count <- ifelse(grepl(" 1 ", results_df_filtered$scenario), "1", "50000")

# Reorder the scenarios by power (calculate average power for ordering)
scenario_order <- results_df_filtered %>%
  group_by(scenario) %>%
  summarize(avg_power = mean(power)) %>%
  arrange(desc(avg_power)) %>%
  pull(scenario)

results_df_filtered$scenario <- factor(results_df_filtered$scenario, levels = scenario_order)

# Reorder PVE factor to control facet order: 2%, 10%, 5%
results_df_filtered$pve <- factor(results_df_filtered$pve, levels = c(0.02, 0.1, 0.05))

# Create the plot
ggplot(results_df_filtered, aes(x = sample_size, y = power, color = scenario, 
                                shape = phenotype_count, linetype = line_type, 
                                group = scenario)) +
  
  # Add vertical reference lines at 500, 1K, 2K
  geom_vline(xintercept = c(500, 1000, 2000), linetype = "dotted", 
             color = "grey30", linewidth = 0.6) +
  
  # Add 80% power reference line
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  
  # Add lines and points with shapes
  geom_line(linewidth = 1.5) +
  geom_point(size = 3.5) +
  
  # Create 2x2 grid layout (10% will be below 2%)
  facet_wrap(~ pve, labeller = labeller(pve = pve_labels), ncol = 2) +
  
  # Labels (no title)
  labs(
    x = "Sample Size (N)",
    y = "Statistical Power",
    color = "Study Design",
    shape = "Number of Phenotypes",
    linetype = "Population Type"
  ) +
  
  # Use log scale for x-axis
  scale_x_log10(
    limits = c(100, 5000),
    breaks = sample_sizes,
    labels = scales::comma
  ) +
  
  # Set y-axis limits and expand to start at 0
  scale_y_continuous(limits = c(0, 1.01), expand = c(0, 0)) +
  
  # Use your custom colors
  scale_color_manual(values = my_colors) +
  
  # Set line types: dashed for 1, solid for 50000
  scale_linetype_manual(values = c("MPP" = "dotdash", "Outbred" = "solid")) +
  
  # Set shapes: triangle (17) for MPP, circle (16) for Outbred
  scale_shape_manual(values = c("1" = 17, "50000" = 16)) +
  
  theme_classic(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.box = "vertical",
    legend.box.just = "left",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  
  # Combine shape and linetype into one guide row
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, order = 1),
    shape = guide_legend(order = 2, title.position = "left"),
    linetype = guide_legend(order = 2, title.position = "left")
  )

# Export 16x12

