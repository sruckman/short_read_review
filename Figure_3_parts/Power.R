# --- 1. Setup: Define Coverage and Power Calculations ---

# Create a sequence of sequencing coverage values
coverage <- seq(1, 105, by = 1)

# Calculate power for the first scenario (SNP-by-SNP)
# This curve rises more slowly.
power_scenario1 <- 1 - exp(-0.05 * coverage)

# Calculate power for a second, more powerful scenario (e.g., Haplotype-based)
# This curve rises more quickly to illustrate a more efficient analysis.
power_scenario2 <- 1 - exp(-0.15 * coverage)


# --- 2. Create the Plot ---

# First, plot the base line for Scenario 1
plot(coverage, power_scenario1,
     type = 'l',
     col = 'dodgerblue',
     lwd = 2,
     xlab = "Sequencing Coverage",
     ylab = "Statistical Power",
     ylim = c(0, 1.0)) # Set y-axis from 0 to 1

# Next, add the second line for Scenario 2 to the existing plot
lines(coverage, power_scenario2,
      col = 'orchid3',
      lwd = 2)


# --- 3. Add Helpful Elements for Interpretation ---

# Add a horizontal line at 80% power, a common benchmark
abline(h = 0.8, col = "grey80", lty = 2)

# Add a legend to explain which line is which
legend("bottomright",
       legend = c("Scenario 1: SNP-by-SNP", "Scenario 2: Haplotype-Based"),
       col = c("dodgerblue", "orchid3"),
       lwd = 2,
       bty = "n") # bty = "n" removes the box around the legend

#Export 6x6