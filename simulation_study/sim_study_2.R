################################################################################
################################# SIMULATION 2 #################################
################################################################################
# Simulation study for joint processes evaluating cross-repetition effects.
# It compares parameter estimates between separable models (fitted on individual
# networks separately) and joint models (fitted on combined networks).
################################################################################

# load required libraries
library(mgcv)
library(dplyr)
library(ggplot2)
library(reshape2)

# load helper functions
source("functions_simu.R")

# set random seed for reproducibility
set.seed(666)

# set parameters
n_simu           <- 1000       # number of simulations
n_students       <- 8          # number of nodes
tot_n_events     <- 200        # total number of events
lambda_0         <- 0.3        # baseline hazard
beta_true        <- 0.5        # true cross-repetition coefficient
delta_rep        <- 0.5        # time decay
n_non_events     <- 1          # number of non-events

# construct the risk set (all directed dyads excluding self-loops)
students <- 1:n_students
dyads    <- expand.grid(sender = students, receiver = students)
dyads    <- dyads[dyads$sender != dyads$receiver, ]
n        <- nrow(dyads)

est_results <- data.frame()



# ------------------------------ SIMULATION LOOP ----------------------------- #

for(sim in 1:n_simu) {
  cat("Simulation", sim, "of", n_simu, "\n")
  
  # initialize event histories for both networks
  history_1 <- list(time=numeric(), sender=integer(), receiver=integer())
  history_2 <- list(time=numeric(), sender=integer(), receiver=integer())
  
  events_1 <- list()
  events_2 <- list()
  
  t <- 0
  
  # generate events until the combined count reaches the target total (tot_n_events) 
  while((length(events_1) + length(events_2)) < tot_n_events) {
    
    # calculate cross-repetition covariate for each network
    cross_1 <- compute_cross_rep(history_2, t, delta_rep)
    cross_2 <- compute_cross_rep(history_1, t, delta_rep)
    
    # calculate hazard rates per dyad
    hazards_1 <- lambda_0 * exp(beta_true * cross_1)
    hazards_2 <- lambda_0 * exp(beta_true * cross_2)
    
    total_hazard_1 <- sum(hazards_1)
    total_hazard_2 <- sum(hazards_2)
    total_system_hazard <- total_hazard_1 + total_hazard_2
    
    if(total_system_hazard == 0) break
    
    # advance time using an exponential waiting time distribution
    dt <- rexp(1, rate = total_system_hazard)
    t  <- t + dt
    
    # determine which network generates the event
    p_net1 <- total_hazard_1 / total_system_hazard
    
    if(runif(1) < p_net1) {
      # event on Network 1
      dyad_idx <- sample(1:n, size = 1, prob = hazards_1 / total_hazard_1)
      s <- dyads$sender[dyad_idx]
      r <- dyads$receiver[dyad_idx]
      
      history_1$time     <- c(history_1$time, t)
      history_1$sender   <- c(history_1$sender, s)
      history_1$receiver <- c(history_1$receiver, r)
      
      events_1[[length(events_1) + 1]] <- list(
        dyad_idx    = dyad_idx,
        time        = t,
        cross_event = cross_1[dyad_idx],
        cross_all   = cross_1
      )
    } else {
      # event on Network 2
      dyad_idx <- sample(1:n, size = 1, prob = hazards_2 / total_hazard_2)
      s <- dyads$sender[dyad_idx]
      r <- dyads$receiver[dyad_idx]
      
      history_2$time     <- c(history_2$time, t)
      history_2$sender   <- c(history_2$sender, s)
      history_2$receiver <- c(history_2$receiver, r)
      
      events_2[[length(events_2) + 1]] <- list(
        dyad_idx    = dyad_idx,
        time        = t,
        cross_event = cross_2[dyad_idx],
        cross_all   = cross_2
      )
    }
  }
  
  # format event pairs for case-control logistic regression
  df_1 <- make_pairs_cross(events_1)
  df_2 <- make_pairs_cross(events_2)
  
  if(!is.null(df_1) && !is.null(df_2)) {
    df_joint <- rbind(df_1, df_2)
    
    # fit separable and joint models using logistic GAMs
    fit_sep_1 <- gam(y ~ -1 + delta_cross, family = binomial, data = df_1)
    fit_sep_2 <- gam(y ~ -1 + delta_cross, family = binomial, data = df_2)
    fit_joint <- gam(y ~ -1 + delta_cross, family = binomial, data = df_joint)
    
    # store coefficient estimates and standard errors
    est_results <- rbind(est_results, data.frame(
      sim        = sim,
      beta_sep_1 = coef(fit_sep_1)["delta_cross"],
      beta_sep_2 = coef(fit_sep_2)["delta_cross"],
      beta_joint = coef(fit_joint)["delta_cross"],
      se_sep_1   = summary(fit_sep_1)$se["delta_cross"],
      se_sep_2   = summary(fit_sep_2)$se["delta_cross"],
      se_joint   = summary(fit_joint)$se["delta_cross"]
    ))
  }
}





# ------------------------- Plots & Visualization ---------------------------- #

# plot with beta estimates
beta_df <- melt(est_results, id.vars = "sim", 
                measure.vars = c("beta_sep_1", "beta_sep_2", "beta_joint"),
                variable.name = "Model", value.name = "Beta")

plot_beta <- ggplot(beta_df, aes(x = Model, y = Beta)) +
  geom_boxplot(fill = "grey", width = 0.6, outlier.shape = NA) +
  geom_hline(yintercept = beta_true, color = "red", linetype = "dashed", linewidth = 0.8) +
  scale_x_discrete(labels = c(
    beta_sep_1 = "separable (net 1)",
    beta_sep_2 = "separable (net 2)",
    beta_joint = "joint (net 1 + net 2)" 
  )) +
  labs(
    y = expression(beta ~ "estimate"),
    x = ""
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.margin = margin(10, 10, 10, 10)
  ) + coord_cartesian(ylim = c(-0.25, 1.25))



# plot with standard errors
se_df <- melt(est_results, id.vars = "sim", 
              measure.vars = c("se_sep_1", "se_sep_2", "se_joint"),
              variable.name = "Model", value.name = "SE")

plot_se <- ggplot(se_df, aes(x = Model, y = SE)) +
  geom_boxplot(fill = "grey", width = 0.5, outlier.shape = NA) +
  scale_x_discrete(labels = c(
    se_sep_1 = "Separable (Net 1)",
    se_sep_2 = "Separable (Net 2)",
    se_joint = "Joint (Net 1 + Net 2)"
  )) +
  labs(
    y = "Standard error",  
    x = "",
    title = ""
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  ) + coord_cartesian(ylim = c(-0.35, 1.2))


print(plot_beta)
#print(plot_se)

#ggsave("plot_beta_simu2.pdf", plot = plot_beta, width = 10, height = 7.5, units = "in", device = cairo_pdf)



# sanity check: joint variance
b1 <- beta_df[1:n_simu, 3]
b2 <- beta_df[(n_simu + 1):(2 * n_simu), 3]
b3 <- beta_df[(2 * n_simu + 1):(3 * n_simu), 3]

v1 <- var(b1)
v2 <- var(b2)
v3 <- var(b3)

print(v1)
print(v2)
print(v3)

v_joint <- (1/v1 + 1/v2)^-1
print(v_joint)