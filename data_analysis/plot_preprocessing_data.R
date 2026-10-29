################################################################################
##################### DESCRIPTIVE ANALYSIS OF GOSSIP DATA ###################### 
################################################################################

library(dplyr)
library(ggplot2)
library(pbapply)
library(ggrepel) 

# source preprocessing script and helper functions
source("data_analysis\\pre_processing_data.R")
source("data_analysis\\functions_data.R")

df <- all_risksets

# map raw 'time' month values (1, 5, 14, 26) to sequential waves (1, 2, 3, 4)
original_times <- sort(unique(df$time))
time_map <- setNames(1:length(original_times), original_times)
df$wave <- time_map[as.character(df$time)]
all_risksets$wave <- df$wave


# -----------------------------------------------------
# GOSSIPING FREQUENCY PER MONTH ACROSS NET1 and NET2
# -----------------------------------------------------
# Calculate monthly gossip rates by filtering active events (z == 1) and grouping 
# by network type and wave

freq_wave_month <- df %>%
  filter(z == 1) %>%
  group_by(network, wave) %>%
  summarise(events = n(),
    months = first(duration_months),
    .groups = "drop"
  ) %>%
  mutate(network = factor(
      network,
      levels = c(1, 2),
      labels = c("self-reported gossip", "perceived gossip")
    )
  ) %>%
  mutate(rate_per_month = events / months
  )


p_freq <- ggplot( freq_wave_month,
  aes(x = factor(wave), y = rate_per_month, fill = network)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6, color = "black") +
  scale_fill_manual(values = c("self-reported gossip" = "grey70", "perceived gossip" = "grey40")
  ) +
  labs(x = "wave", y = "gossiping frequency per month", fill = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
    limits = c(0, max(freq_wave_month$rate_per_month, na.rm = TRUE) * 1.1)
  )+
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    legend.position = "top",
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_freq)

# save plot
#ggsave("gossip_real_perc_monthly.pdf", plot = p_freq, width = 9, height = 7)




# ------------------------------------------------------------------------------
# SCATTER PLOT: GOSSIPING FREQUENCY PER CLASS 
# ------------------------------------------------------------------------------
# load student metadata dataset
students_df <- read.csv(
  "data_analysis/gossipdata/RECENS_Wired_into_each_other-Data_files/RECENS_students_dataset_w1234.csv"
)

# lookup tables (wave -> idcode -> class)
students_lookup <- list(
  "1" = setNames(students_df$class_1, students_df$idcode),
  "2" = setNames(students_df$class_2, students_df$idcode),
  "3" = setNames(students_df$class_3, students_df$idcode),
  "4" = setNames(students_df$class_4, students_df$idcode)
)

# map ClassGroup across all risk sets
all_risksets$ClassGroup <- pbmapply(
  get_group_class_fast,
  all_risksets$sender,
  all_risksets$receiver,
  all_risksets$wave
)
all_risksets$ClassGroup <- factor(all_risksets$ClassGroup)

# summarize events per network and class
class_summary_net <- all_risksets %>%
  filter(z == 1) %>%                               
  group_by(network, ClassGroup) %>%
  summarise(n_events = n(), .groups = "drop")

# a scatter plot to compare self-reported and perceived gossip frequency by class
class_scatter <- class_summary_net %>%
  mutate(
    network = recode(
      network,
      `1` = "self-reported gossip",
      `2` = "perceived gossip"
    )
  ) %>%
  tidyr::pivot_wider(
    names_from = network,
    values_from = n_events
  )

p_class_scatter <- ggplot(
  class_scatter,
  aes(
    x = `self-reported gossip`,
    y = `perceived gossip`
  )
) +
  geom_point(
    shape = 16,
    size = 1.8,
    fill = "black",
    color = "black"
  )+
  geom_text_repel(
    aes(label = ClassGroup),
    size = 3.2,
    box.padding = 0.15,       
    point.padding = 0.08,     
    max.overlaps = Inf,      
    segment.color = "grey50",
    segment.size = 0.3,
    fontface = "bold"
  ) +
  labs(
    x = "self-reported gossip frequency",
    y = "perceived gossip frequency"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.05))
  ) +
  theme_classic(base_size = 20) +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_class_scatter)
# save plot
#ggsave("gossip_real_perc_classes.pdf", plot = p_class_scatter, width = 9, height = 7)




# -------------------------------------------
# GENDER DYAD ANALYSIS (F-F, F-M, M-F, M-M)
# -------------------------------------------
# Combine network covariate datasets and categorize ties based on sender-receiver gender pairs
gender_pairs_summary <- bind_rows(
  riskset_network1_cov %>% mutate(network = "self-reported gossip"),
  riskset_network2_cov %>% mutate(network = "perceived gossip")
) %>% 
  filter(z == 1) %>%
  mutate(network = factor(
      network, levels = c("self-reported gossip", "perceived gossip")
    ),
    gender_pair = case_when(
      girl_ego == 1 & girl_alter == 1 ~ "F-F",
      girl_ego == 0 & girl_alter == 0 ~ "M-M",
      girl_ego == 1 & girl_alter == 0 ~ "F-M",
      girl_ego == 0 & girl_alter == 1 ~ "M-F"
    ),
    gender_pair = factor(
      gender_pair,
      levels = c("F-F","M-M","F-M","M-F")
    )
  ) %>%
  group_by(network, gender_pair) %>%
  summarise(events = n(), .groups = "drop")

gender_pairs_summary

# create a bar chart comparing event counts across gender dyad categories
p_gender_pairs <- ggplot(gender_pairs_summary, aes(x = gender_pair, y = events, fill = network)
) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6, color = "black") +
  scale_fill_manual(values = c("self-reported gossip" = "grey70", "perceived gossip" = "grey40")
  ) +
  labs(y = "gossiping frequency", x="", fill = NULL) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(
      size = 20,
      margin = margin(r = 10)),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "top",
    legend.text = element_text(size = 20),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

p_gender_pairs
# save plot
#ggsave("gossip_real_perc_gender.pdf", plot = p_gender_pairs, width = 9, height = 7)
