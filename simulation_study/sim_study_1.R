################################################################################
################################# SIMULATION 1 #################################
################################################################################
# Simulation study for a joint process generating two simultaneous networks.
# It compares models with between-type cross-reciprocity versus within-type 
# repetition using logistic regression (GAMs). Evaluates performance via AIC/BIC.
################################################################################

# load required libraries
library(mgcv)
library(mgcv)
library(dplyr)
library(ggplot2)
library(reshape2)

# load helper functions
source("functions_simu.R")

# set random seed for reproducibility
set.seed(666) 

# set parameters
n_simu       <- 50      # number of simulation runs
tot_n_events <- 2000    # total combined events across both networks
n_students   <- 8       # number of nodes (students)

# baseline hazard rates for the two processes
lambda_0_1    <- 0.15 
lambda_0_2    <- 0.17 

# true beta coefficients for cross-reciprocity
beta_cross_1  <- 0.9 
beta_cross_2  <- 0.7

# temporal decay parameters for reciprocity and repetition
delta_rec    <- 1 
delta_rep    <- delta_rec/2
n_non_events <- 1


# construct the risk set (all directed dyads excluding self-loops)
students <- 1:n_students
dyads <- expand.grid(sender = students, receiver = students)
dyads <- dyads[dyads$sender != dyads$receiver, ]
n <- nrow(dyads)

# initialize data frame for results
results <- data.frame()
beta_estimates <- data.frame()
fits <- list(
  cross = vector("list", n_simu),
  rep   = vector("list", n_simu),
  both  = vector("list", n_simu)
)



# ------------------------------ SIMULATION LOOP ----------------------------- #

for(sim in 1:n_simu){
  
  cat("\nSimulation", sim, "of", n_simu, "\n")
  
  # initialize event histories for both networks
  history_1 <- list(time=numeric(), sender=integer(), receiver=integer())
  history_2 <- list(time=numeric(), sender=integer(), receiver=integer())
  
  cross_1 <- rep(0, n)
  cross_2 <- rep(0, n)
  
  t <- 0
  events_1 <- list()
  events_2 <- list()
  
  # generate events until the combined count reaches the target total (tot_n_events)
  while((length(events_1) + length(events_2)) < tot_n_events){
    
    # compute current hazards based on history/covariates
    hazard_1 <- lambda_0_1 * exp(beta_cross_1 * cross_1)
    hazard_2 <- lambda_0_2 * exp(beta_cross_2 * cross_2)
    
    total_intensity <- sum(hazard_1) + sum(hazard_2)
    if(total_intensity == 0) break
    
    # advance time using an exponential waiting time distribution
    t <- t + rexp(1, total_intensity)
    
    # determine which network generates the event
    if(runif(1) < sum(hazard_1)/total_intensity){
      
      # event occurs in Network 1
      dyad_idx <- sample(1:n, 1, prob = hazard_1)
      event <- dyads[dyad_idx, ]
      
      cross_current <- compute_cross(history_2, t, delta_rec)
      rep_current   <- compute_rep_decay(history_1, t, delta_rep)
      
      events_1[[length(events_1)+1]] <- list(
        dyad_idx = dyad_idx,
        time = t,
        cross_event = cross_current[dyad_idx],
        cross_all = cross_current,
        rep_event = rep_current[dyad_idx],
        rep_all = rep_current
      )
      
      history_1$time     <- c(history_1$time, t)
      history_1$sender   <- c(history_1$sender, event$sender)
      history_1$receiver <- c(history_1$receiver, event$receiver)
      cross_1 <- cross_current
      
    } else {
      
      # event occurs in Network 2
      dyad_idx <- sample(1:n, 1, prob = hazard_2)
      event <- dyads[dyad_idx, ]
      
      cross_current <- compute_cross(history_1, t, delta_rec)
      rep_current   <- compute_rep_decay(history_2, t, delta_rep)
      
      events_2[[length(events_2)+1]] <- list(
        dyad_idx = dyad_idx,
        time = t,
        cross_event = cross_current[dyad_idx],
        cross_all = cross_current,
        rep_event = rep_current[dyad_idx],
        rep_all = rep_current
      )
      
      history_2$time     <- c(history_2$time, t)
      history_2$sender   <- c(history_2$sender, event$sender)
      history_2$receiver <- c(history_2$receiver, event$receiver)
      cross_2 <- cross_current
    }
  }
  # format event pairs for statistical modelling
  pairs_1 <- make_pairs(events_1)
  pairs_2 <- make_pairs(events_2)
  
  all_pairs <- rbind(
    cbind(pairs_1, network = 1),
    cbind(pairs_2, network = 2)
  )
  
  # fit candidate models using logistic GAMs (binomial family)
  # 1. model with only cross-reciprocity
  fit_cross_1 <- gam(y ~ -1 + delta_cross, family = binomial, data = pairs_1)
  fit_cross_2 <- gam(y ~ -1 + delta_cross, family = binomial, data = pairs_2)
  
  # 2. model with only within-repetition
  fit_rep_1 <- gam(y ~ -1 + delta_rep, family = binomial, data = pairs_1)
  fit_rep_2 <- gam(y ~ -1 + delta_rep, family = binomial, data = pairs_2)
  
  # 3. model with both covariates
  fit_both_1 <- gam(y ~ -1 + delta_cross + delta_rep, family = binomial, data = pairs_1)
  fit_both_2 <- gam(y ~ -1 + delta_cross + delta_rep, family = binomial, data = pairs_2)
  
  # save information criteria (AIC and BIC) for comparison
  results <- rbind(results, 
                   data.frame(
                     sim = sim,
                     model = c("cross","rep","both"),
                     AIC = c(AIC(fit_cross_1), AIC(fit_rep_1), AIC(fit_both_1)),
                     BIC = c(BIC(fit_cross_1), BIC(fit_rep_1), BIC(fit_both_1))
                   ))
  
  # save beta coefficient estimates 
  beta_estimates <- rbind(beta_estimates,
                          data.frame(
                            sim = sim,
                            beta_cross_only = coef(fit_cross_1)["delta_cross"],
                            beta_cross_both = coef(fit_both_1)["delta_cross"]
                          ))
  
  # store model fit
  fits$cross[[sim]] <- fit_cross_1
  fits$rep[[sim]]   <- fit_rep_1
  fits$both[[sim]]  <- fit_both_1
  
}





# ------------------------- Plots & Visualization ---------------------------- #

# prepare model factor levels
results$model <- factor(results$model, levels = c("cross", "rep", "both"))

# Note: For the paper, we focus strictly on comparing the model with 
# cross-reciprocity vs. the model with within-repetition (dropping 'both').
res <- subset(results, model %in% c("cross", "rep"))
res$model <- factor(res$model, levels = c("cross", "rep"))


# BIC plot
plot_BIC <- ggplot(res, aes(x = model, y = BIC, fill = model)) +
  geom_boxplot(width = 0.6) +
  theme_minimal() +
  labs(
    y = "BIC",
    x = "",
    title = ""
  ) +
  scale_x_discrete(labels = c(
    cross = "Between-type reciprocity",
    rep   = "Within-type repetition"
  )) +
  scale_fill_manual(values = c("grey", "grey")) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(size = 20),   
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20)  
  )


# AIC plot
y_limits <- res %>%
  group_by(model) %>%
  summarise(
    baffo_min = boxplot.stats(AIC)$stats[1],
    baffo_max = boxplot.stats(AIC)$stats[5]
  ) %>%
  summarise(
    min_tot = min(baffo_min),
    max_tot = max(baffo_max)
  ) %>%
  as.numeric()

plot_AIC <- ggplot(res, aes(x = model, y = AIC, fill = model)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  coord_cartesian(ylim = y_limits) +
  #theme_minimal() +
  labs(y = "AIC", x = "") +
  scale_x_discrete(labels = c(
    cross = "between-type reciprocity",
    rep   = "within-type repetition"
  )) +
  scale_fill_manual(values = c("grey", "grey")) +
  theme_classic(base_size = 20)+
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.margin = margin(10, 10, 10, 10)
  )


# print and save the AIC plot
print(plot_AIC)
#ggsave("plot_AIC_simu1.pdf", plot = plot_AIC, width = 10, height = 7.5, units = "in", device = cairo_pdf)
