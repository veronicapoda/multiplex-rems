################################################################################
############################# DATA PREPROCESSING ###############################
################################################################################
# Data preprocessing for the multiplex relational event model.
# Processes raw survey waves, constructs self-reported (Net 1) and perceived (Net 2)
# gossip networks, builds risk sets, and integrates exogenous student covariates.


library(dplyr)

source("utils\\covariates.R")
source("data_analysis\\functions_data.R")



# --------------------
# DATA AND WAVE SETUP
# --------------------
waves <- c("Wave1", "Wave2", "Wave3", "Wave4")
wave_map <- setNames(1:length(waves), waves)
folders <- c("1000", "2000", "3000", "4000", "5000", "6000", "7000")

data_path <- paste0("data_analysis/gossipdata/RECENS_Wired_into_each_other-Data_files/")
network_path <- paste0(data_path, "RECENS_networks_dataset_w1234/")
students_path <- paste0(data_path, "RECENS_students_dataset_w1234.csv")


# -------------------------------------------------------------
# CONSTRUCTION OF SELF-REPORTED AND PERCEIVED GOSSIP NETWORKS
# -------------------------------------------------------------
# initialize data frame for all observed gossip events across layers
all_gossip_events <- data.frame(
  sender = character(),
  receiver = character(),
  time = integer(),
  network = integer(),
  stringsAsFactors = FALSE
)

# loop over survey waves, school folders, and classes to extract network matrices
for (wave in waves) {
  
  wave_num <- wave_map[[wave]]
  
  cat("\n\nWave:", wave, "\n")
  
  codes <- get_network_codes(wave)
  
  for (folder in folders) {
    
    folder_path <- file.path(network_path, wave, folder)
    
    if (!dir.exists(folder_path)) {
      next
    }
    
    classes_in_folder <- list.files(path = folder_path, full.names = FALSE)
    
    for (class in classes_in_folder) {
      
      # exclude classes starting with 9
      if (grepl("^9", class)) {
        next
      }
      
      class_path <- file.path(folder_path, class)
      
      if (length(list.files(class_path)) == 0) {
        next
      }
      
      # self-reported gossip (Network 1)
      gossip_self <- process_matrix(
        folder_path,
        class,
        codes["self"],
        network_id = 1
      )
      
      # perceived gossip (Network 2)
      gossip_perceived <- process_matrix(
        folder_path,
        class,
        codes["perceived"],
        network_id = 2
      )
      
      # assign wave
      if (nrow(gossip_self) > 0) {
        gossip_self$time <- wave_num
      }
      
      if (nrow(gossip_perceived) > 0) {
        gossip_perceived$time <- wave_num
      }
      
      # combine both layers for the current class
      gossip_all <- bind_rows(
        gossip_self,
        gossip_perceived
      )
      
      # remove self-loops
      gossip_all <- gossip_all %>%
        filter(sender != receiver)
      
      if (nrow(gossip_all) > 0) {
        all_gossip_events <- bind_rows(
          all_gossip_events,
          gossip_all
        )
      }
    }
  }
}


cat(
  "\nTotal observed gossip events:",
  nrow(all_gossip_events),
  "\n"
)


# ---------------------------------
# SEPARATE THE TWO NETWORK LAYERS
# ---------------------------------

network1_events <- all_gossip_events %>%
  filter(network == 1)

network2_events <- all_gossip_events %>%
  filter(network == 2)

cat(
  "Self-reported gossip events:",
  nrow(network1_events),
  "\n"
)

cat(
  "Perceived gossip events:",
  nrow(network2_events),
  "\n"
)



# -------------------------
# TEMPORAL REPRESENTATION
# -------------------------
# duration of each survey wave in months
durations <- c(2, 6, 12, 12)

# cumulative wave starting times
start_time <- cumsum(c(0, durations[-length(durations)]))

# midpoints of observation intervals
mid_time <- start_time + durations / 2

# wave midpoints:
# wave 1 = 1, wave 2 = 5, wave 3 = 14, wave 4 = 26


# -----------------------
# RISK-SET CONSTRUCTION
# -----------------------
all_risksets <- data.frame()

for (wave in waves) {
  
  wave_num <- wave_map[[wave]]
  
  cat("\nWave:", wave, "\n")
  
  codes <- get_network_codes(wave)
  
  for (folder in folders) {
    
    folder_path <- file.path(
      network_path,
      wave,
      folder
    )
    
    if (!dir.exists(folder_path)) {
      next
    }
    
    classes_in_folder <- list.files(folder_path)
    
    for (class in classes_in_folder) {
      
      if (grepl("^9", class)) {
        next
      }
      
      class_path <- file.path(
        folder_path,
        class
      )
      
      if (length(list.files(class_path)) == 0) {
        next
      }
      
      # read a reference matrix to retrieve all active students in the class
      self_path <- file.path(
        folder_path,
        class,
        substr(codes["self"], 1, 2),
        paste0(class, "_", codes["self"], ".csv")
      )
      
      if (!file.exists(self_path)) {
        next
      }
      
      mat <- as.matrix(
        read.csv(
          self_path,
          row.names = 1
        )
      )
      all_ids <- rownames(mat)
      
      # filter observed events for the current wave
      events_self <- network1_events %>%
        filter(time == wave_num)
      
      events_perceived <- network2_events %>%
        filter(time == wave_num)
      
      # generate risk sets for both network layers
      rs1 <- generate_risk_set(
        all_ids,
        wave_num,
        network_id = 1,
        events_self
      )
      rs2 <- generate_risk_set(
        all_ids,
        wave_num,
        network_id = 2,
        events_perceived
      )
      
      all_risksets <- bind_rows(
        all_risksets,
        rs1,
        rs2
      )
    }
  }
}


cat(
  "\nTotal risk-set observations:",
  nrow(all_risksets),
  "\n"
)



# ------------------------------
# SEPARATE RISK SETS BY NETWORK
# ------------------------------
riskset_network1 <- all_risksets %>%
  filter(network == 1)

riskset_network2 <- all_risksets %>%
  filter(network == 2)


# ----------------------
# EXOGENOUS COVARIATES
# ----------------------
students_df <- read.csv(
  students_path
)

# gender
students_df$gender_final <- apply(
  students_df[, c(
    "gender_1",
    "gender_2",
    "gender_3",
    "gender_4"
  )],
  1,
  function(x) {
    x <- x[!is.na(x)]
    
    if (length(x) == 0) {
      return(NA)
    }
    
    x[1]
  }
)

# recode gender: 1 = male -> 0, 2 = female -> 1
students_df$gender_final <- ifelse(
  students_df$gender_final == 2,
  1,
  ifelse(
    students_df$gender_final == 1,
    0,
    NA
  )
)

# create gender lookup vector
gender_vec <- students_df$gender_final
names(gender_vec) <- students_df$idcode

# map gender covariates onto risk sets
riskset_network1$girl_ego <-
  gender_sender(
    riskset_network1,
    gender_vec
  )

riskset_network1$girl_alter <-
  gender_receiver(
    riskset_network1,
    gender_vec
  )

riskset_network2$girl_ego <-
  gender_sender(
    riskset_network2,
    gender_vec
  )

riskset_network2$girl_alter <-
  gender_receiver(
    riskset_network2,
    gender_vec
  )


# same-gender indicator
riskset_network1 <- riskset_network1 %>%
  mutate(
    same_gender =
      as.numeric(girl_ego == girl_alter)
  )

riskset_network2 <- riskset_network2 %>%
  mutate(
    same_gender =
      as.numeric(girl_ego == girl_alter)
  )


# filter out observations with missing gender information
riskset_network1 <- riskset_network1 %>%
  filter(
    !is.na(girl_ego),
    !is.na(girl_alter)
  )
riskset_network2 <- riskset_network2 %>%
  filter(
    !is.na(girl_ego),
    !is.na(girl_alter)
  )


# ---------------
# SANITY CHECKS
# ---------------
cat("\n--- Number of students ---\n")
cat(
  "Network 1:",
  length(unique(riskset_network1$sender)),
  "\n"
)
cat(
  "Network 2:",
  length(unique(riskset_network2$sender)),
  "\n"
)

cat("\n--- Number of observed events ---\n")
print(
  table(riskset_network1$z)
)
print(
  table(riskset_network2$z)
)

cat("\n--- Same-gender proportions ---\n")
print(
  prop.table(
    table(riskset_network1$same_gender)
  )
)
print(
  prop.table(
    table(riskset_network2$same_gender)
  )
)



# -----------------------------
# FINAL PREPROCESSED DATASETS
# -----------------------------
riskset_network1_cov <- riskset_network1 %>%
  select(
    sender,
    receiver,
    wave,
    time,
    duration_months,
    offset_log,
    z,
    girl_ego,
    girl_alter,
    same_gender
  )

riskset_network2_cov <- riskset_network2 %>%
  select(
    sender,
    receiver,
    wave,
    time,
    duration_months,
    offset_log,
    z,
    girl_ego,
    girl_alter,
    same_gender
  )
