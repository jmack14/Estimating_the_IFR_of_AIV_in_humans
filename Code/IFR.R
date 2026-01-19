# ============================================================
# Estimating the infection fatality ratio for humans infected with avian influenza viruses
# ============================================================

# ============================================================
# Setup
# ============================================================

setwd("C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/figures")
load("Interpandemic_period_results.RData")

library(patchwork)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)
library(ggrepel)
library(grid)
library(lhs)
library(sensitivity)
library(triangle)

# ============================================================
# pi_i function
# ============================================================

model <- function(parameters) {
  a     <- parameters[1]
  Rstar <- parameters[2]
  R0    <- parameters[3]
  psi   <- parameters[4]
  
  a1 <- (R0 - 1 - a)
  a2 <- (R0 + 1 + a)^2
  P  <- (1 - (1 / Rstar))
  a3 <- 1 + (a * (1 - P))
  
  (a1 + sqrt(a2 - 4 * R0 * a3)) / (2 * R0)
}

# ============================================================
# Run a single scenario
# ============================================================

run_scenario <- function(param_ranges, num_samples) {
  
  # Latin Hypercube Sampling
  lhs_samples <- randomLHS(num_samples, length(param_ranges))
  
  transformed_samples <- data.frame(
    a = qunif(lhs_samples[,1],
              min = param_ranges$a[1],
              max = param_ranges$a[2]),
    Rstar = qtriangle(lhs_samples[,2],
                      a = param_ranges$Rstar[1],
                      b = param_ranges$Rstar[2],
                      c = param_ranges$Rstar[3]),
    R0 = qexp(lhs_samples[,3],
              rate = param_ranges$R0$rate),
    psi = qgamma(lhs_samples[,4],
                 shape = param_ranges$psi$shape,
                 scale = param_ranges$psi$scale)
  )
  
  # Calculate pi_i
  transformed_samples$pi <- apply(transformed_samples, 1, model)
  bar_pi <- mean(transformed_samples$pi)
  
  # Solve for mean number of zoonotic spillovers
  Lambda <- function(mean_nz, Lambda_target) {
    nz <- seq(0, 1e5)
    PNz_nz <- dpois(nz, mean_nz)
    1 - sum(PNz_nz * (1 - bar_pi)^nz) - Lambda_target
  }
  
  Lambda_target_vec <- 1 / transformed_samples$psi
  mean_nz_vec <- numeric(length(Lambda_target_vec))
  
  for (i in seq_along(Lambda_target_vec)) {
    mean_nz_vec[i] <- uniroot(Lambda, interval = c(0, 1e6), Lambda_target = Lambda_target_vec[i])$root
  }
  
  transformed_samples$mean_nz <- mean_nz_vec
  
  # Mean number of human infections
  bar_m <- mean(1 / (1 - transformed_samples$R0))
  transformed_samples$meannh <- bar_m * transformed_samples$mean_nz
  
  # Sort and calculate cumulative probabilities
  cum.prob <- seq(1 / num_samples, 1, 1 / num_samples)
  sort.output <- data.frame(
    meannh = sort(transformed_samples$meannh),
    cum.prob = cum.prob
  )
  
  # Quantiles
  i025 <- min(which(sort.output$cum.prob >= 0.025))
  i25  <- min(which(sort.output$cum.prob >= 0.25))
  i50  <- min(which(sort.output$cum.prob >= 0.5))
  i75  <- min(which(sort.output$cum.prob >= 0.75))
  i975 <- min(which(sort.output$cum.prob >= 0.975))
  
  quarts <- data.frame(
    cprob = c(sort.output$cum.prob[i025], sort.output$cum.prob[i25],
              sort.output$cum.prob[i50], sort.output$cum.prob[i75],
              sort.output$cum.prob[i975]),
    n = c(sort.output$meannh[i025], sort.output$meannh[i25],
          sort.output$meannh[i50], sort.output$meannh[i75],
          sort.output$meannh[i975])
  )
  
  list(
    transformed_samples = transformed_samples,
    sort_output = sort.output,
    quarts = quarts
  )
}

# ============================================================
# Run scenarios
# ============================================================

num_samples <- 1000

res_1 <- run_scenario(
  list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta/theta0, scale = theta0)
  ),
  num_samples
)

res_2 <- run_scenario(
  list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta0/theta, scale = theta)
  ),
  num_samples
)

res_3 <- run_scenario(
  list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta0/theta, scale = 18.5/(delta0/theta))
  ),
  num_samples
)

# ============================================================
# Save results
# ============================================================

write.csv(res_1$transformed_samples, "transformed_samples_38.csv", row.names = FALSE, quote = FALSE)
write.csv(res_2$transformed_samples, "transformed_samples_53.5.csv", row.names = FALSE, quote = FALSE)
write.csv(res_3$transformed_samples, "transformed_samples_18.5.csv", row.names = FALSE, quote = FALSE)

write.csv(res_1$sort_output, "mean_nh_38.csv", row.names = FALSE, quote = FALSE)
write.csv(res_2$sort_output, "mean_nh_53.5.csv", row.names = FALSE, quote = FALSE)
write.csv(res_3$sort_output, "mean_nh_18.5.csv", row.names = FALSE, quote = FALSE)

# ============================================================
# Plot the mean annual number of human infections with AIV
# ============================================================

Figure_2A <- ggplot() +
  geom_line(data = res_1$sort_output, aes(x = log(meannh,10), y = cum.prob, colour = "orange"), size = 1.1) +
  geom_line(data = res_2$sort_output, aes(x = log(meannh,10), y = cum.prob, colour = "blue"), size = 1.1) +
  geom_line(data = res_3$sort_output, aes(x = log(meannh,10), y = cum.prob, colour = "red"), size = 1.1) +
  geom_point(data = res_1$quarts, aes(x = log(n,10), y = cprob)) +
  geom_point(data = res_2$quarts, aes(x = log(n,10), y = cprob)) +
  geom_point(data = res_3$quarts, aes(x = log(n,10), y = cprob)) +
  geom_text_repel(data = res_1$quarts, aes(x = log(n,10), y = cprob, label = round(n,0)), size = 6, hjust = 0, vjust = 1) +
  geom_text_repel(data = res_2$quarts, aes(x = log(n,10), y = cprob, label = round(n,0)), size = 6, hjust = 0, vjust = 1) +
  geom_text_repel(data = res_3$quarts, aes(x = log(n,10), y = cprob, label = round(n,0)), size = 6, hjust = 1, vjust = 0) +
  ylab("Quantiles") +xlab(expression(bold("Mean number of human infections, ") * bolditalic(bar(n))[h]))+
  scale_x_continuous(breaks = seq(3,5), labels = c("1,000","10,000","100,000"), limits = c(3,5)) +
  scale_color_manual(values = c("blue","red","orange"), labels = c("53.5","38","18.5")) +
  labs(color = "Mean interpandemic\nperiod (years)") +
  labs(tag = expression(bold("A")))+
  theme_bw() +
  theme(axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = c(0.8,0.4),
        legend.key.size = unit(1.5,"cm"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0,1))

# ============================================================
# Prepare IFR distributions
# ============================================================

deaths <- 20.6

make_ifr_df <- function(res){
  df <- res$sort_output %>%
    mutate(IFR = deaths / meannh * 10000)
  
  data.frame(
    IFR = sort(df$IFR),
    cum.prob = seq(1 / nrow(df), 1, 1 / nrow(df))
  )
}

ifr_1 <- make_ifr_df(res_1)
ifr_2 <- make_ifr_df(res_2)
ifr_3 <- make_ifr_df(res_3)

# ----------------------------
# Quantiles 
# ----------------------------

get_quarts <- function(ifr_df){
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  idx <- sapply(probs, function(p)
    min(which(ifr_df$cum.prob >= p))
  )
  
  data.frame(
    cprob = ifr_df$cum.prob[idx],
    IFR   = ifr_df$IFR[idx]
  )
}

q1 <- get_quarts(ifr_1)
q2 <- get_quarts(ifr_2)
q3 <- get_quarts(ifr_3)

# ============================================================
# Plot the infection fatality ratios: Figure 2B
# ============================================================

Figure_2B <- ggplot() +
  geom_line(data = ifr_1, aes(x = log(IFR,10), y = cum.prob, colour = "orange"), size = 1.1) +
  geom_line(data = ifr_2, aes(x = log(IFR,10), y = cum.prob, colour = "blue"),   size = 1.1) +
  geom_line(data = ifr_3, aes(x = log(IFR,10), y = cum.prob, colour = "red"),    size = 1.1) +
  
  geom_point(data = q1, aes(x = log(IFR,10), y = cprob)) +
  geom_point(data = q2, aes(x = log(IFR,10), y = cprob)) +
  geom_point(data = q3, aes(x = log(IFR,10), y = cprob)) +
  
  geom_text_repel(data = q1, aes(x = log(IFR,10), y = cprob, label = round(IFR,1)),
                  size = 6, hjust = 0, vjust = 1) +
  geom_text_repel(data = q2, aes(x = log(IFR,10), y = cprob, label = round(IFR,1)),
                  size = 6, hjust = 0, vjust = 1) +
  geom_text_repel(data = q3, aes(x = log(IFR,10), y = cprob, label = round(IFR,1)),
                  size = 6, hjust = 1, vjust = 0) +
  
  ylab("Quantiles") +
  xlab("Infection fatality ratio per 10,000 infections") +
  
  scale_x_continuous(
    breaks = seq(0, 3),
    labels = c("1","10","100","1000"),
    limits = c(0.5, 2.25)
  ) +
  
  scale_color_manual(
    values = c("blue","red","orange"),
    labels = c("53.5","38","18.5")
  ) +
  
  labs(tag = expression(bold("B"))) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",   # ✅ legend removed
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0,1)
  )

# ============================================================
# Combine Figure 2A and 2B 
# ============================================================
Figure_2 <- Figure_2A + Figure_2B +
  plot_layout(ncol = 1, heights = c(1, 1)) 

ggsave("Figure_2.png", width = 12, height = 14, dpi = 300)

# ============================================================
# Plot parameter distributions: Figure A1
# ============================================================

plot_distributions <- function(transformed_samples) {
  options(scipen = 999)
  
  plot_a <- ggplot(transformed_samples, aes(x = a)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(x = "Probability", y = "Count") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                       axis.text = element_text(size = 16, face = "bold"),
                       axis.title = element_text(size = 16, face = "bold")) +
    ggtitle("Probability that a strain of HPAI with pandemic potential emerges, *a*") +
    theme(plot.title = element_markdown())
  
  plot_Rstar <- ggplot(transformed_samples, aes(x = Rstar)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(x = "Probability", y = "Count") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                       axis.text = element_text(size = 16, face = "bold"),
                       axis.title = element_text(size = 16, face = "bold")) +
    ggtitle("Reproduction number of a reassortment HPAI virus in humans, *R**") +
    theme(plot.title = element_markdown())
  
  plot_R0 <- ggplot(transformed_samples, aes(x = R0)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(x = "Probability", y = "Count") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                       axis.text = element_text(size = 16, face = "bold"),
                       axis.title = element_text(size = 16, face = "bold")) +
    ggtitle("Reproduction number of a HPAI virus in humans prior to evolutionary change, *R*<sub>0") +
    theme(plot.title = element_markdown())
  
  return((plot_a)/(plot_Rstar)/(plot_R0))
}

Figure_A1 <- plot_distributions(res_1$transformed_samples)

ggsave("Figure_A1.png", Figure_A1, height = 8, width = 10)