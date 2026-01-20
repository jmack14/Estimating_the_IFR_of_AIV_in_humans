# ============================================================
# Risk of human infection with HPAI from bird migration and risk of coinfection with human influenza
# ============================================================

# ============================================================
# Load libraries
# ============================================================

library(ggplot2)
library(patchwork)

# ============================================================
# Seasonal risk 
# ============================================================

season_data <- data.frame(
  Month = factor(
    c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
    levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  ),
  MigrationRisk = factor(
    c("Low","Med","High","High","Med","Low","Low","Low","Med","High","Med","Low"),
    levels = c("Low","Med","High")
  ),
  CoinfectionRisk = factor(
    c("High","High","Med","Med","Low","Low","Low","Low","Low","Med","Med","High"),
    levels = c("Low","Med","High")
  )
)

# ============================================================
# Plotting helper
# ============================================================

risk_plot <- function(data, risk_var, title, fill_color) {
  ggplot(data, aes(x = Month, y = !!risk_var)) +
    geom_bar(stat = "identity", fill = fill_color) +
    labs(y = "Risk", x = NULL, title = title) +
    theme_bw() +
    theme(
      axis.text  = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 16, face = "bold")
    )
}

# ============================================================
# Generate plots
# ============================================================

Migration <- risk_plot(
  season_data,
  quote(MigrationRisk),
  title = "Migration",
  fill_color = "blue"
)

HI <- risk_plot(
  season_data,
  quote(CoinfectionRisk),
  title = "Coinfection",
  fill_color = "red"
)

# ============================================================
# Combine plots: Figure 5
# ============================================================

Figure_5 <- Migration / HI
Figure_5

ggsave("Figure_5.png", Figure_5, height = 6, width = 8)
