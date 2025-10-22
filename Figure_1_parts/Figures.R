#Figures for Review Paper

library(ggplot2)
library(scales)
library(dplyr)
library(ggrepel)
library(tidyr)

df <- cost_per_mb[c(1:78),]
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")



ggplot(df, aes(x = Date, y = `Cost per Mb`)) +
  geom_line(color = "dodgerblue", linewidth = 1.2) +
  scale_y_log10(
    limits = c(0.001, 10000),
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000),
    labels = c("$0.001", "$0.01", "$0.10", "$1", "$10", "$100", "$1,000", "$10,000")
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0.05, 0)
  ) +
  labs(
    x = "Year",
    y = "Cost per Megabase (log scale)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

#Export 6x7

ggplot(df, aes(x = Date, y = `Cost per Read`)) +
  geom_line(color = "dodgerblue", linewidth = 1.2) +
  scale_y_log10(
    limits = c(0.0000001, 1),
    breaks = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1),
    labels = c("$0.0000001", "$0.000001", "$0.00001", "$0.0001","$0.001", "$0.01", "$0.10", "$1")
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0.05, 0)
  ) +
  labs(
    x = "Year",
    y = "Cost per Read (log scale)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )

#Export 6x7

df_long <- df %>% pivot_longer( cols = c(`Cost per Mb`, `Cost per Read`), names_to = "Metric", values_to = "Value" )
ggplot(df_long, aes(x = Date, y = Value, color = Metric)) + 
  geom_line(linewidth = 1.2) + facet_wrap(~ Metric, scales = "free_y", ncol = 1) + 
  scale_y_log10( labels = scales::dollar_format() ) + 
  scale_x_date( date_breaks = "2 years", date_labels = "%Y" ) + 
  scale_color_manual( values = c("Cost per Mb" = "dodgerblue", "Cost per Read" = "orchid") ) + 
  labs( x = "Year", y = "Cost (log scale)") + 
  theme_minimal(base_size = 14) + 
  theme( axis.text.x = element_text(angle = 45, hjust = 1), 
         panel.grid = element_blank(), 
         axis.line = element_line(color = "black", linewidth = 0.8), 
         axis.ticks = element_line(color = "black", linewidth = 0.8), 
         strip.text = element_text(face = "bold", size = 14), legend.position = "none" )
