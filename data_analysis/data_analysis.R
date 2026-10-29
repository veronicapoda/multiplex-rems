################################################################################
##################### DYNAMIC MULTIPLEX NETWORK ANALYSIS #######################
################################################################################
# Empirical analysis for the multiplex relational event model on gossip data.
# Computes endogenous covariates across network layers (self-reported & perceived).
# Estimates linear, smooth and random effects using generalized additive mixed models.

library(mgcv)
library(dplyr)
library(mgcViz)
library(ggplot2)
library(pbapply)
library(xtable)

source("utils/covariates.R")
source("data_analysis/functions_data.R")
source("data_analysis/pre_processing_data.R")


# ------------------------------------------------------------------------------
# Network 1 - self-reported gossip
# ------------------------------------------------------------------------------

# compute endogenous covariates for Network 1 (averaged history type)
endogenous_covs_1 <- c(
  "reciprocity",
  "transitive.closure", 
  "cyclic.closure", 
  "sender.balance", 
  "receiver.balance"
)
covs_network1 <- compute_endogenous_covs(riskset_network1, endogenous_covs_1, type = "avg")

# compute specific past-history covariates for Network 1
cov_repetition <- compute_endogenous_covs(riskset_network1, "repetition", type = "past")
cov_sender_out_degree <- compute_endogenous_covs(riskset_network1,"sender.out.degree", type = "past")
cov_receiver_in_degree <- compute_endogenous_covs(riskset_network1,"receiver.in.degree", type = "past")

covs_network1 <- cbind(
  cov_repetition,
  cov_sender_out_degree,
  cov_receiver_in_degree,
  covs_network1
)

# note on perceived gossip network:
# - raw data: "sender -> receiver" means the sender thinks the receiver gossips about them.
# - our model: "sender -> receiver" means the receiver thinks the sender gossips about them (flipped perspective).
# because of this flip, we must interpret directional effects accordingly.

endogenous_covs_2 <- c(
  "repetition",
  "reciprocity", 
  "sender.out.degree",
  "receiver.in.degree",
  "sender.in.degree",
  "receiver.out.degree")
covs_network2 <- compute_endogenous_covs(riskset_network2, endogenous_covs_2, type = "avg")



# merge covariates into respective risk set
riskset_network1_cov <- bind_cols(
  riskset_network1 %>% dplyr::select(sender, receiver, wave, time, offset_log, z, girl_ego, girl_alter, same_gender),
  covs_network1
)
riskset_network2_cov <- riskset_network2 %>%
  dplyr::select(sender, receiver, wave, time, offset_log, z, girl_ego, girl_alter, same_gender) %>%
  bind_cols(covs_network2)

# construct combined data for model 1
data_model1 <- riskset_network1_cov %>%
  dplyr::select(sender, receiver, wave, time, offset_log, girl_ego, girl_alter, same_gender, z, everything()) %>%
  rename_with(~ paste0(.x, "1"),
              -c(sender, receiver, wave, time, offset_log, girl_ego, girl_alter, same_gender, z)) %>%
  mutate(
    z1 = z,
    girl_ego = girl_ego,
    girl_alter = girl_alter,
    same_gender = same_gender
  ) %>%
  dplyr::select(-z) %>%
  bind_cols(
    riskset_network2_cov %>%
      dplyr::select(-sender, -receiver, -wave, -time, -offset_log, -girl_ego, -girl_alter, -same_gender) %>%
      rename_with(~ paste0(.x, "2"))
  ) %>%
  dplyr::select(sender, receiver, wave, time, offset_log, girl_ego, girl_alter, same_gender, z1, everything())

# re-orient directional covariates for the perceived network layer in model 1
data_model1 <- data_model1 %>%
  mutate(
    old_sender_out2   = sender.out.degree2,
    old_receiver_in2  = receiver.in.degree2,
    old_sender_in2    = sender.in.degree2,
    old_receiver_out2 = receiver.out.degree2
  ) %>%
  mutate(
    sender.in.degree2   = old_receiver_out2,
    receiver.in.degree2 = old_sender_out2,  
    sender.out.degree2    = old_receiver_in2,  
    receiver.out.degree2  = old_sender_in2    
  ) %>%
  dplyr::select(-old_sender_out2, -old_receiver_in2, -old_sender_in2, -old_receiver_out2)


reciprocity2_past <- compute_endogenous_covs(riskset_network2, "repetition", type = "past")
data_model1$delayed_retaliation2 <- reciprocity2_past$repetition

repetition2 <- compute_endogenous_covs(riskset_network2, "reciprocity", type = "avg")
data_model1$repetition2 <- repetition2$reciprocity

data_model1$immediate_retaliation2 <- data_model1$z2
data_model1$delayed_repetition1 <- data_model1$repetition1

data_model1$repetition1 <- NULL
data_model1$z2 <- NULL
data_model1$reciprocity2 <- NULL
data_model1$repetition2 <- NULL
data_model1$past_sender.out.degree1 <- data_model1$sender.out.degree1
data_model1$sender.out.degree1 <- NULL
data_model1$past_receiver.in.degree1 <- data_model1$receiver.in.degree1
data_model1$receiver.in.degree1 <- NULL


# apply square root transformation to degree variables
data_model1 <- data_model1 %>%
  mutate(
    past_sender.out.degree1 = sqrt(past_sender.out.degree1),
    past_receiver.in.degree1 = sqrt(past_receiver.in.degree1),
    sender.out.degree2 = sqrt(sender.out.degree2),
    receiver.in.degree2 = sqrt(receiver.in.degree2),
    sender.in.degree2 = sqrt(sender.in.degree2),
    receiver.out.degree2 = sqrt(receiver.out.degree2)
  )

# define smooth and linear terms for model 1 formula
smooth_covariates <- c(
  "time", 
  "delayed_repetition1",
  "reciprocity1",
  "past_sender.out.degree1", 
  "past_receiver.in.degree1",
  "transitive.closure1", 
  "cyclic.closure1",
  "sender.balance1",
  "receiver.balance1",
  "delayed_retaliation2",
  "receiver.out.degree2",
  "sender.in.degree2"
)
smooth_terms <- sapply(smooth_covariates, gen_smooth, data = data_model1)
linear_terms <- c("girl_ego", 
                  "girl_alter",
                  "same_gender",
                  "immediate_retaliation2")

formula_full <- as.formula(
  paste0(
    "z1 ~ ",            
    paste(c(linear_terms, smooth_terms), collapse=" + "),
    " + offset(offset_log)"
  )
)

# initial fit without random effects
fit1_smooth <- bam(
  formula_full,
  family = binomial(link = "cloglog"),
  data = data_model1,
  select = TRUE,
  method = "fREML",
  discrete = TRUE
)

summary(fit1_smooth)
plot(fit1_smooth, pages = 1, scale = 0, rug = TRUE)
AIC(fit1_smooth)



# -----------------------------
# Model 1 with random effects
# -----------------------------
# account for unobserved heterogeneity at class, sender, and receiver levels

students_df <- read.csv("data_analysis/gossipdata/RECENS_Wired_into_each_other-Data_files/RECENS_students_dataset_w1234.csv")
students_lookup <- list(
  "1" = setNames(students_df$class_1, students_df$idcode),
  "2" = setNames(students_df$class_2, students_df$idcode),
  "3" = setNames(students_df$class_3, students_df$idcode),
  "4" = setNames(students_df$class_4, students_df$idcode)
)
original_times <- sort(unique(data_model1$time))  
time_map <- setNames(1:length(original_times), original_times)
data_model1$wave <- time_map[as.character(data_model1$time)]
unique(data_model1$wave)

# class level
data_model1$ClassGroup <- pbmapply(
  get_group_class_fast,
  data_model1$sender,
  data_model1$receiver,
  data_model1$wave
)
data_model1$ClassGroup <- factor(data_model1$ClassGroup)
table(data_model1$ClassGroup, useNA = "ifany")

# receiver level
data_model1$ReceiverGroup <- factor(data_model1$receiver)

# sender level
data_model1$SenderGroup <- factor(data_model1$sender)

formula_full_re <- update(formula_full, . ~ . + 
                            s(ClassGroup, bs = "re") + 
                            s(ReceiverGroup, bs = "re") + 
                            s(SenderGroup, bs = "re"))

# re-fit with random effects 
fit1_smooth_re <- bam(
  formula_full_re,
  family = binomial(link = "cloglog"),
  data = data_model1,
  select = TRUE,
  method = "fREML",
  discrete = TRUE
)

summary(fit1_smooth_re)
plot(fit1_smooth_re, pages = 1, scale = 0, rug = TRUE)
AIC(fit1_smooth_re)
gam.vcomp(fit1_smooth_re)

plot(fit1_smooth_re, select = 4, rug = TRUE, scale=0)
plot(fit1_smooth_re, select = 4, rug = TRUE, jit = TRUE, scale=0)

# plot select smooth effects for model 1 using mgcViz
p_model1 <- getViz(fit1_smooth_re)

dev.off()
#pdf("sender.out.degree_S.pdf", width = 8, height = 6)
p1 <- plot(sm(p_model1, 4)) +
  l_fitLine() +
  l_ciLine() +
  l_rug() +
  theme_classic(base_size = 20) +
  labs(x = "sender out-degree (S)", y = "log-hazard contribution (S)") +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),   
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),   
    axis.text.x = element_text(size = 20),    
    axis.text.y = element_text(size = 20),   
    plot.margin  = margin(10, 10, 10, 10)
  )
print(p1)
dev.off()



#pdf("retaliation_P.pdf", width = 8, height = 6)
p2 <- plot(sm(p_model1, 10)) +
  l_fitLine() +
  l_ciLine() +
  l_rug() +
  theme_classic(base_size = 20) +
  labs(x = "delayed retaliation (P)", y = "log-hazard contribution (S)") +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),  
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),   
    axis.text.x = element_text(size = 20),  
    axis.text.y = element_text(size = 20),   
    plot.margin  = margin(10, 10, 10, 10)
  )
print(p2)
dev.off()








# ------------------------------------------------------------------------------
# Network 2 - perceived gossip
# ------------------------------------------------------------------------------

endogenous_covs_1 <- c(
  "sender.out.degree",
  "receiver.in.degree", 
  "sender.in.degree",
  "receiver.out.degree",
  "transitive.closure",
  "cyclic.closure", 
  "sender.balance",
  "receiver.balance"
  )
covs_network1 <- compute_endogenous_covs(riskset_network1, endogenous_covs_1, type = "avg") 

covs_rep <- compute_endogenous_covs(
  riskset_network2, 
  c("repetition",
    "sender.out.degree",
    "receiver.in.degree"
    ), 
  type = "past"
  )
covs_other <- compute_endogenous_covs(riskset_network2,"reciprocity", type = "avg")
covs_network2 <- cbind(covs_rep, covs_other)

riskset_network1_cov <- NULL
riskset_network2_cov <- NULL
riskset_network1_cov <- bind_cols(
  riskset_network1 %>% dplyr::select(sender, receiver, wave, time, offset_log, z, girl_ego, girl_alter, same_gender),
  covs_network1
)
riskset_network2_cov <- riskset_network2 %>%
  dplyr::select(sender, receiver, wave, time, offset_log, z, girl_ego, girl_alter, same_gender) %>%
  bind_cols(covs_network2)

# construct dataset for model 2
data_model2 <- riskset_network2_cov %>%
  dplyr::select(sender, receiver, wave, time, offset_log, girl_ego, girl_alter, same_gender, z, dplyr::everything()) %>%
  rename_with(~ paste0(.x, "2"), -c(sender, receiver, wave, time, offset_log, girl_ego, girl_alter, same_gender, z)) %>%
  mutate(z2 = z) %>%
  dplyr::select(-z) %>%
  bind_cols(
    riskset_network1_cov %>%
      dplyr::select(-sender, -receiver, -wave, -time, -offset_log, -girl_ego, -girl_alter, -same_gender) %>%
      rename_with(~ paste0(.x, "1"))
  ) %>%
  dplyr::select(sender, receiver, wave, time, offset_log, girl_ego, girl_alter, same_gender, z2, dplyr::everything())

repetition1_past <- compute_endogenous_covs(riskset_network1, "repetition", type = "past") 
repetition1_past <- repetition1_past$repetition

data_model2$delayed_retaliation1 <- repetition1_past
data_model2$delayed_repetition2 <- data_model2$repetition2
data_model2$repetition2 <- NULL

# re-map directional network 1 covariates for model 2
data_model2$new_receiver.out.degree1 <- data_model2$sender.out.degree1
data_model2$new_sender.in.degree1    <- data_model2$receiver.in.degree1
data_model2$new_receiver.in.degree1  <- data_model2$sender.in.degree1
data_model2$new_sender.out.degree1   <- data_model2$receiver.out.degree1

data_model2$receiver.out.degree1 <- data_model2$new_receiver.out.degree1
data_model2$sender.in.degree1    <- data_model2$new_sender.in.degree1
data_model2$receiver.in.degree1  <- data_model2$new_receiver.in.degree1
data_model2$sender.out.degree1   <- data_model2$new_sender.out.degree1

data_model2$new_receiver.out.degree1 <- NULL
data_model2$new_sender.in.degree1    <- NULL
data_model2$new_receiver.in.degree1  <- NULL
data_model2$new_sender.out.degree1   <- NULL

data_model2$new_receiver.out.degree2 <- data_model2$sender.out.degree2
data_model2$new_sender.in.degree2    <- data_model2$receiver.in.degree2
data_model2$past_receiver.out.degree2 <- data_model2$new_receiver.out.degree2
data_model2$past_sender.in.degree2    <- data_model2$new_sender.in.degree2

data_model2$new_receiver.out.degree2 <- NULL
data_model2$new_sender.in.degree2    <- NULL
data_model2$sender.out.degree2       <- NULL
data_model2$receiver.in.degree2      <- NULL

# define covariates modelled through non-linear smooth effects
smooth_covariates_2 <- c(
  "time",
  "delayed_repetition2", 
  "reciprocity2",
  "delayed_retaliation1",
  "receiver.in.degree1", 
  "sender.out.degree1",
  "receiver.out.degree1",
  "sender.in.degree1", 
  "cyclic.closure1", 
  "transitive.closure1", 
  "sender.balance1",
  "receiver.balance1"
)
data_model2 <- data_model2 %>% rename(immediate_retaliation1 = z1)

# define covariates included as linear effects
linear_terms_2 <- c("girl_ego", "girl_alter", "same_gender", "immediate_retaliation1")

# apply square root transformation to degree variables
smooth_terms_2 <- sapply(smooth_covariates_2, gen_smooth_sqrt, data = data_model2)

formula_full_2 <- as.formula(
  paste0("z2 ~ ",
         paste(c(linear_terms_2, smooth_terms_2), collapse=" + "),
         " + offset(offset_log)")
)
formula_full_2


# fit model 2 without random effects
fit2_smooth <- bam(
  formula_full_2,
  family = binomial(link = "cloglog"),
  data = data_model2,
  select = TRUE,
  method = "fREML",    
  discrete = TRUE
)

summary(fit2_smooth)
plot(fit2_smooth, pages = 1, scale = 0, rug = TRUE)
AIC(fit2_smooth)


# -----------------------------
# Model 2 with random effects
# -----------------------------
# account for unobserved heterogeneity at class, sender, and receiver levels

students_df <- read.csv("data_analysis/gossipdata/RECENS_Wired_into_each_other-Data_files/RECENS_students_dataset_w1234.csv")
students_lookup <- list(
  "1" = setNames(students_df$class_1, students_df$idcode),
  "2" = setNames(students_df$class_2, students_df$idcode),
  "3" = setNames(students_df$class_3, students_df$idcode),
  "4" = setNames(students_df$class_4, students_df$idcode)
)
original_times <- sort(unique(data_model2$time))  
time_map <- setNames(1:length(original_times), original_times)
data_model2$wave <- time_map[as.character(data_model2$time)]
unique(data_model2$wave)

# class level
data_model2$ClassGroup <- pbmapply(
  get_group_class_fast,
  data_model2$sender,
  data_model2$receiver,
  data_model2$wave
)
data_model2$ClassGroup <- factor(data_model2$ClassGroup)
table(data_model2$ClassGroup, useNA = "ifany")

# receiver level
data_model2$ReceiverGroup <- factor(data_model2$receiver)

# sender level
data_model2$SenderGroup <- factor(data_model2$sender)

formula_full2_re <- update(formula_full_2, . ~ . + 
                             s(ClassGroup, bs = "re") + 
                             s(ReceiverGroup, bs = "re") + 
                             s(SenderGroup, bs = "re"))

# fit model 2 with random effects
fit2_smooth_re <- bam(
  formula_full2_re,
  family = binomial(link = "cloglog"),
  data = data_model2,
  select = TRUE,
  method = "fREML",
  discrete = TRUE
)

# data analysis & plot
summary(fit2_smooth_re)
plot(fit2_smooth_re, pages = 1, scale = 0, rug = TRUE)
AIC(fit2_smooth_re)
gam.vcomp(fit2_smooth_re)

# plot smooth effects for model 2 using mgcViz
p_model2 <- getViz(fit2_smooth_re)
dev.off()
#pdf("reciprocity_P.pdf", width = 8, height = 6)
p1 <- plot(sm(p_model2, 3)) +
  l_fitLine() +
  l_ciLine() +
  l_rug() +
  theme_classic(base_size = 20) +
  labs(x = "reciprocity (P)", y = "log-hazard contribution (P)") +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),  
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),   
    axis.text.x = element_text(size = 20),    
    axis.text.y = element_text(size = 20),   
    plot.margin  = margin(10, 10, 10, 10)
  )
print(p1)
dev.off()

#pdf("receiver.in.degree_S.pdf", width = 8, height = 6)
p2 <- plot(sm(p_model2, 5)) +
  l_fitLine() +
  l_ciLine() +
  l_rug() +
  theme_classic(base_size = 20) +
  labs(x = "receiver in-degree (S)", y = "log-hazard contribution (P)") +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),  
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),   
    axis.text.x = element_text(size = 20),   
    axis.text.y = element_text(size = 20),  
    plot.margin  = margin(10, 10, 10, 10)
  )
print(p2)
dev.off()







# -----------------------
# LATEX EXPORT UTILITIES
# -----------------------
# convert p-values into standard significance star notation
sig_stars <- function(p) {
  as.character(symnum(p, corr = FALSE, na = FALSE,
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", "·", "")))
}

# format p-values: use LaTeX scientific notation for extremely small values
format_p <- function(p) {
  ifelse(p < 2e-16, "$<2\\times10^{-16}$", formatC(p, format = "e", digits = 2))
}

# export GAM fixed and smooth tables to LaTeX format
export_gam_latex <- function(fit) {
  s_fit <- if (inherits(fit, "summary.gam")) fit else summary(fit)
  
  # fixed effects
  fixed <- as.data.frame(s_fit$p.table[, 1:4])
  colnames(fixed) <- c("Estimate", "Std. Error", "z value", "p.raw")
  fixed$p <- format_p(fixed$p.raw)
  fixed$stars <- sig_stars(fixed$p.raw)
  rownames(fixed) <- gsub("_", " ", rownames(fixed))
  fixed <- fixed[, c("Estimate", "Std. Error", "z value", "p", "stars")]
  colnames(fixed)[4] <- "p-value"
  cat("% --- FIXED EFFECTS ---\n")
  print(xtable(fixed, digits = 3), include.rownames = TRUE, only.contents = TRUE,
        sanitize.text.function = identity, type = "latex")
  
  # smooth effects 
  smooth <- as.data.frame(s_fit$s.table[, 1:4])
  colnames(smooth) <- c("edf", "Ref.df", "Chi.sq", "p.raw")
  smooth$Ref.df <- ifelse(is.na(smooth$Ref.df), "", formatC(smooth$Ref.df, format = "f", digits = 0))
  smooth$edf    <- formatC(smooth$edf, format = "f", digits = 3)
  smooth$Chi.sq <- formatC(smooth$Chi.sq, format = "f", digits = 3)
  smooth$p      <- format_p(smooth$p.raw)
  smooth$stars  <- sig_stars(smooth$p.raw)
  rownames(smooth) <- gsub("_", " ", rownames(smooth))
  smooth <- smooth[, c("edf", "Ref.df", "Chi.sq", "p", "stars")]
  colnames(smooth)[4] <- "p-value"
  cat("% --- SMOOTH EFFECTS ---\n")
  print(xtable(smooth, digits = 3), include.rownames = TRUE, only.contents = TRUE,
        sanitize.text.function = identity, type = "latex")
}

export_gam_latex(fit1_smooth_re)
export_gam_latex(fit2_smooth_re)