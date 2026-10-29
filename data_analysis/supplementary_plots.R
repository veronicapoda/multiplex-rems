################################################################################
########################### SUPPLEMENTARY PLOTS ################################
################################################################################

library(mgcViz)
library(ggplot2)

source("data_analysis/data_analysis.R")



# ------------------------------------------------------------------------------
# Model 1 - self-reported gossip network
# ------------------------------------------------------------------------------
p_model1 <- getViz(fit1_smooth_re)

smooth_labels_1 <- c(
  "time" = "time",
  "delayed_repetition1" = "repetition (S)",
  "reciprocity1" = "reciprocity (S)",
  "past_sender.out.degree1" = "sender out-degree (S)",
  "past_receiver.in.degree1" = "receiver in-degree (S)",
  "transitive.closure1" = "transitive closure (S)",
  "cyclic.closure1" = "cyclic closure (S)",
  "sender.balance1" = "sender balance (S)",
  "receiver.balance1" = "receiver balance (S)",
  "delayed_retaliation2" = "delayed retaliation (P)",
  "receiver.out.degree2" = "receiver out-degree (P)",
  "sender.in.degree2" = "sender in-degree (P)"
)

supp_indices_1 <- setdiff(
  seq_along(smooth_covariates),
  c(4, 10)
)

plot_list_1 <- lapply(supp_indices_1, function(i) {
  var <- smooth_covariates[i]
  label_x <- if (!is.na(smooth_labels_1[var])) smooth_labels_1[var] else var
  
  plot(sm(p_model1, i)) +
    l_fitLine() +
    l_ciLine() +
    l_rug() +
    theme_classic(base_size = 20) +
    labs(
      x = label_x,
      y = "log-hazard contribution (S)"
    ) +
    theme(
      axis.title.x = element_text(size = 20, margin = margin(t = 10)),
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),
      axis.text.x  = element_text(size = 20),
      axis.text.y  = element_text(size = 20),
      plot.margin  = margin(30, 50, 30, 50)
    )
})

pdf(
  "supplementary_smooths_model1.pdf",
  width = 16,
  height = 22
)
do.call(gridPrint, c(plot_list_1, ncol = 2))
dev.off()
graphics.off()




# ------------------------------------------------------------------------------
# Model 2 - perceived gossip network
# ------------------------------------------------------------------------------
p_model2 <- getViz(fit2_smooth_re)

smooth_labels_2 <- c(
  "time" = "time",
  "delayed_repetition2" = "repetition (P)",
  "reciprocity2" = "reciprocity (P)",
  "delayed_retaliation1" = "delayed retaliation (S)",
  "receiver.in.degree1" = "receiver in-degree (S)",
  "sender.out.degree1" = "sender out-degree (S)",
  "receiver.out.degree1" = "receiver out-degree (S)",
  "sender.in.degree1" = "sender in-degree (S)",
  "cyclic.closure1" = "cyclic closure (S)",
  "transitive.closure1" = "transitive closure (S)",
  "sender.balance1" = "sender balance (S)",
  "receiver.balance1" = "receiver balance (S)"
)

supp_indices_2 <- setdiff(
  seq_along(smooth_covariates_2),
  c(3, 5)
)

plot_list_2 <- lapply(supp_indices_2, function(i) {
  
  var <- smooth_covariates_2[i]
  label_x <- if (!is.na(smooth_labels_2[var])) smooth_labels_2[var] else var
  
  plot(sm(p_model2, i)) +
    l_fitLine() +
    l_ciLine() +
    l_rug() +
    theme_classic(base_size = 20) +
    labs(
      x = label_x,
      y = "log-hazard contribution (P)"
    ) +
    theme(
      axis.title.x = element_text(size = 20, margin = margin(t = 10)),
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),
      axis.text.x  = element_text(size = 20),
      axis.text.y  = element_text(size = 20),
      plot.margin  = margin(30, 50, 30, 50)
    )
})

pdf(
  "supplementary_smooths_model2.pdf",
  width = 16,
  height = 22
)
do.call(gridPrint, c(plot_list_2, ncol = 2))
dev.off()
