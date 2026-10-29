################################################################################
################################ UTILITY FUNCTIONS #############################
################################################################################

# return specific file codes for self-reported and perceived networks per wave
get_network_codes <- function(wave) {
  
  if (wave == "Wave1") {
    return(c(self = "33_3", perceived = "33_4"))
  }
  
  if (wave == "Wave2") {
    return(c(self = "32_3_2h", perceived = "32_4_2h"))
  }
  
  if (wave == "Wave3") {
    return(c(self = "39_3_3h", perceived = "39_4_3h"))
  }
  
  if (wave == "Wave4") {
    return(c(self = "38_3_4h", perceived = "38_4_4h"))
  }
  
  stop("Unknown wave.")
}


# read a class adjacency matrix and convert it into an event data frame
process_matrix <- function(
    folder_path,
    class,
    code,
    network_id
) {
  
  code_prefix <- substr(code, 1, 2)
  
  path <- file.path(
    folder_path,
    class,
    code_prefix,
    paste0(class, "_", code, ".csv")
  )
  
  if (!file.exists(path)) {
    return(
      data.frame(
        sender = character(),
        receiver = character(),
        network = integer(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  mat <- as.matrix(
    read.csv(path, row.names = 1)
  )
  
  rows <- rownames(mat)
  cols <- colnames(mat)
  
  idx <- which(mat == 1, arr.ind = TRUE)
  
  if (nrow(idx) == 0) {
    return(
      data.frame(
        sender = character(),
        receiver = character(),
        network = integer(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  data.frame(
    sender = rows[idx[, 1]],
    receiver = sub("^X", "", cols[idx[, 2]]),
    network = network_id,
    stringsAsFactors = FALSE
  )
}


# generate complete risk set with temporal offsets and event indicators
generate_risk_set <- function(
    ids,
    wave_num,
    network_id,
    events
) {
  
  rs <- expand.grid(
    sender = ids,
    receiver = ids,
    stringsAsFactors = FALSE
  ) %>%
    filter(sender != receiver)
  
  rs$wave <- wave_num
  
  # midpoint representation of time
  rs$time <- mid_time[wave_num]
  
  # duration of observation interval
  rs$duration_months <- durations[wave_num]
  
  # log duration offset
  rs$offset_log <- log(
    rs$duration_months
  )
  
  # event indicator (1 if dyad experienced an event, 0 otherwise)
  if (nrow(events) > 0) {
    
    event_keys <- paste(
      events$sender,
      events$receiver,
      sep = "_"
    )
    
    risk_keys <- paste(
      rs$sender,
      rs$receiver,
      sep = "_"
    )
    
    rs$z <- as.integer(
      risk_keys %in% event_keys
    )
    
  } else {
    
    rs$z <- 0
    
  }
  
  rs$network <- network_id
  
  rs
}


# generate standard GAM smooth term string with adaptive basis dimension (k)
gen_smooth <- function(var, data) {
  k_val <- min(length(unique(data[[var]])), 10)
  sprintf("s(%s, k=%d)", var, k_val)
}


# return the class of a student given their ID and wave time
get_class_for_student_fast <- function(student_id, time) {
  cls <- students_lookup[[as.character(time)]][as.character(student_id)]
  
  if (is.na(cls) || is.null(cls) || length(cls) == 0) return(NA)
  return(cls)
}


# return the most frequent class among senders and receivers for group-level effects
get_group_class_fast <- function(sender, receiver, time) {
  all_ids <- c(sender, receiver)
  classes <- sapply(all_ids, get_class_for_student_fast, time = time)
  classes <- na.omit(classes)
  if (length(classes) == 0) return(NA)
  class_mode <- as.numeric(names(sort(table(classes), decreasing = TRUE)[1]))
  return(class_mode)
}


# generate GAM smooth term string incorporating a square-root transformation for degree variables
gen_smooth_sqrt <- function(var, data) {
  if (var %in% c("sender.out.degree1", "receiver.in.degree1", "sender.in.degree1", 
                 "receiver.out.degree1", "past_sender.in.degree2", "past_receiver.out.degree2")) {
    new_var <- paste0("sqrt(", var, ")")
  } else {
    new_var <- var
  }
  k_val <- min(length(unique(data[[var]])), 10)
  sprintf("s(%s, k=%d)", new_var, k_val)
}
