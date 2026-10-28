################################################################################ 
################################ COVARIATES ####################################
################################################################################
# Functions for computing exogenous and endogenous covariates used in the 
# Multiplex Relational Event Model. 
################################################################################


# -------------------------- EXOGENOUS COVARIATES ---------------------------- #

# sender gender
gender_sender <- function(RiskSet, gender_vec){
  sapply(RiskSet$sender, function(s) ifelse(!is.na(gender_vec[as.character(s)]),
                                            gender_vec[as.character(s)], NA))
}

# receiver gender 
gender_receiver <- function(RiskSet, gender_vec){
  sapply(RiskSet$receiver, function(r) ifelse(!is.na(gender_vec[as.character(r)]),
                                              gender_vec[as.character(r)], NA))
}

# age 
age <- function(RiskSet) {
  13 + RiskSet$time  # Wave 1 → 14, Wave 2 → 15, etc.
}






# ------------------------ ENDOGENOUS COVARIATES ----------------------------- #

# ---------------------
# monadic effects 
# ---------------------

# sender out-degree: 
# number of events sent by the sender prior to the current event
sender_out_degree <- function(history, event) {
  if (is.null(history) || event$time <= 1) return(0)
  sum(history$sender == event$sender & history$time <= event$time)
}

# receiver in-degree: 
# number of events received by the receiver prior to the current event
receiver_in_degree <- function(history, event) {
  if (is.null(history) || event$time <= 1) return(0)
  sum(history$receiver == event$receiver & history$time <= event$time)
}

# receiver out-degree: 
# number of events sent by the receiver prior to the current event
receiver_out_degree <- function(history, event) {
  if (is.null(history) || event$time <= 1) return(0)
  sum(history$sender == event$receiver & history$time <= event$time)
}

# sender in-degree: 
# number of events received by the sender prior to the current event
sender_in_degree <- function(history, event) {
  if (is.null(history) || event$time <= 1) return(0)
  sum(history$receiver == event$sender & history$time <= event$time)
}




# ------------------ 
# dyadic effects
# ------------------ 

# repetition: 
# number of previous events between the same sender-receiver pair
repetition <- function(history, event) {
  if (is.null(history) || event$time <= 1) return(0)
  sum(history$sender == event$sender & history$receiver == event$receiver &
        history$time <= event$time)
}

# reciprocity: 
# number of previous events in the reverse direction (receiver -> sender)
reciprocity <- function(history, event) {
  if (is.null(history) || event$time <= 1) return(0)
  sum(history$sender == event$receiver & history$receiver == event$sender &
        history$time <= event$time)
}




# ----------------- 
# triadic effects 
# ----------------- 

# transitive closure: 
# counts patterns of the form: sender -> intermediary -> receiver 
transitive_closure <- function(history, event) {
  time <- event$time
  s <- event$sender
  r <- event$receiver
  
  if (is.null(history) || time <= 1) return(0)
  
  tc_sum <- 0
  counted_nodes <- c()
  
  # past events s -> a, where a != s and a != r
  past_idx <- which(history$sender == s & history$receiver != s & history$receiver != r & history$time <= time)
  
  for (i in past_idx) {
    a <- history$receiver[i]
    hy_deg_s_a <- length(unique(history$time[history$sender == s & history$receiver == a & history$time <= time]))
    hy_deg_a_r <- length(unique(history$time[history$sender == a & history$receiver == r & history$time <= time]))
    
    if (!(a %in% counted_nodes)) {
      tc_sum <- tc_sum + min(hy_deg_s_a, hy_deg_a_r)
      counted_nodes <- c(counted_nodes, a)
    }
  }
  
  return(tc_sum)
}

# cyclic closure 
# counts patterns of the form: receiver -> intermediary -> sender 
cyclic_closure <- function(history, event) {
  time <- event$time
  s <- event$sender
  r <- event$receiver
  
  if (is.null(history) || time <= 1) return(0)
  
  cc_sum <- 0
  counted_nodes <- c()
  
  # past events a -> s, where a != s and a != r
  past_idx <- which(history$receiver == s & history$sender != s & history$sender != r & history$time <= time)
  
  for (i in past_idx) {
    a <- history$sender[i]
    hy_deg_a_s <- length(unique(history$time[history$sender == a & history$receiver == s & history$time <= time]))
    hy_deg_r_a <- length(unique(history$time[history$sender == r & history$receiver == a & history$time <= time]))
    
    if (!(a %in% counted_nodes)) {
      cc_sum <- cc_sum + min(hy_deg_a_s, hy_deg_r_a)
      counted_nodes <- c(counted_nodes, a)
    }
  }
  
  return(cc_sum)
}

# sender balance
# counts overlap between incoming ties: a -> sender and a -> receiver
sender_balance <- function(history, event) {
  time <- event$time
  s <- event$sender
  r <- event$receiver
  
  if (is.null(history) || time <= 1) return(0)
  
  sb_sum <- 0
  counted_nodes <- c()
  
  past_idx <- which(history$receiver == s & history$sender != s & history$sender != r & history$time <= time)
  
  for (i in past_idx) {
    a <- history$sender[i]
    hy_deg_a_s <- length(unique(history$time[history$sender == a & history$receiver == s & history$time <= time]))
    hy_deg_a_r <- length(unique(history$time[history$sender == a & history$receiver == r & history$time <= time]))
    
    if (!(a %in% counted_nodes)) {
      sb_sum <- sb_sum + min(hy_deg_a_s, hy_deg_a_r)
      counted_nodes <- c(counted_nodes, a)
    }
  }
  
  return(sb_sum)
}

# receiver balance 
# counts overlap between outgoing ties: sender -> a and receiver -> a
receiver_balance <- function(history, event) {
  time <- event$time
  s <- event$sender
  r <- event$receiver
  
  if (is.null(history) || time <= 1) return(0)
  
  rb_sum <- 0
  counted_nodes <- c()
  
  past_idx <- which(history$sender == s & history$receiver != s & history$receiver != r & history$time <= time)
  
  for (i in past_idx) {
    a <- history$receiver[i]
    hy_deg_s_a <- length(unique(history$time[history$sender == s & history$receiver == a & history$time <= time]))
    hy_deg_r_a <- length(unique(history$time[history$sender == r & history$receiver == a & history$time <= time]))
    
    if (!(a %in% counted_nodes)) {
      rb_sum <- rb_sum + min(hy_deg_s_a, hy_deg_r_a)
      counted_nodes <- c(counted_nodes, a)
    }
  }
  
  return(rb_sum)
}







# --------------------------------
# GENERIC COVARIATE FUNCTION
# -------------------------------- 
# function to compute either exogenous and endogenous covariates 
compute_covariates <- function(RiskSet, 
                               CovariateName, 
                               History = NULL, 
                               Event = NULL, 
                               gender_map = NULL) {
  
  
  # ---------------------- 
  # endogenous covariates 
  # ----------------------
  if (CovariateName %in% c("sender.out.degree",
                           "receiver.in.degree",
                           "sender.in.degree",
                           "receiver.out.degree", 
                           "repetition", 
                           "reciprocity",
                           "transitive.closure",
                           "cyclic.closure",
                           "sender.balance",
                           "receiver.balance")) {
    
    if (is.null(History) || is.null(Event)) stop("Endogenous covariates require History and Event")
    
    if (CovariateName == "sender.out.degree") return(sender_out_degree(History, Event))
    if (CovariateName == "receiver.in.degree") return(receiver_in_degree(History, Event))
    if (CovariateName == "sender.in.degree") return(sender_in_degree(History, Event))
    if (CovariateName == "receiver.out.degree") return(receiver_out_degree(History, Event))
    if (CovariateName == "repetition") return(repetition(History, Event))
    if (CovariateName == "reciprocity") return(reciprocity(History, Event))
    if (CovariateName == "transitive.closure") return(transitive_closure(History, Event))
    if (CovariateName == "cyclic.closure") return(cyclic_closure(History, Event))
    if (CovariateName == "sender.balance") return(sender_balance(History, Event))
    if (CovariateName == "receiver.balance") return(receiver_balance(History, Event))
    
  } else {
    
    # ---------------------- 
    # exogenous covariates 
    # ----------------------
    if (CovariateName == "girl.ego") return(gender_sender(RiskSet, gender_map))
    if (CovariateName == "girl.alter") return(gender_receiver(RiskSet, gender_map))
    if (CovariateName == "age") return(age(RiskSet))
    
    # direct risk set lookup  
    if (!(CovariateName %in% colnames(RiskSet))) stop("Covariate not found in dataset: ", CovariateName)
    return(RiskSet[[CovariateName]])
  }
}






# ----------------------------------------
# ENDOGENOUS COVARIATES EVALUATION TYPES
# ----------------------------------------
# options: 
# - "past" : history strictly before the current event 
# - "curr" : history including the current event 
# - "avg"  : average of past and current values 
compute_endogenous_covs <- function(risk_set_df, 
                                    covariate_names,
                                    type = c("avg", "curr", "past")) {
  
  type <- match.arg(type)
  z_data <- risk_set_df$z
  n <- nrow(risk_set_df)
  
  # history containing observed events only
  history_data <- list(time = numeric(0),
                       receiver = character(0),
                       sender = character(0))
  
  cov_past    <- lapply(covariate_names, function(x) numeric(n))
  cov_current <- lapply(covariate_names, function(x) numeric(n))
  cov_avg     <- lapply(covariate_names, function(x) numeric(n))
  
  names(cov_past)    <- covariate_names
  names(cov_current) <- covariate_names
  names(cov_avg)     <- covariate_names
  
  for (i in 1:n) {
    
    e <- list(
      time     = risk_set_df$time[i],
      receiver = risk_set_df$receiver[i],
      sender   = risk_set_df$sender[i]
    )
    
    # ------------------------------------------------------------ 
    # past history: events strictly before the current event time 
    # ------------------------------------------------------------
    past_idx <- which(history_data$time < e$time)
    history_past <- list(
      time     = history_data$time[past_idx],
      receiver = history_data$receiver[past_idx],
      sender   = history_data$sender[past_idx]
    )
    
    # -------------------------------------------
    # current history: include the current event
    # -------------------------------------------
    if (z_data[i] == 1) {
      history_current <- list(
        time     = c(history_data$time,     e$time),
        receiver = c(history_data$receiver, e$receiver),
        sender   = c(history_data$sender,   e$sender)
      )
    } else {
      history_current <- history_data
    }
     
    # compute covariates 
    for (cov_name in covariate_names) {
      
      # past
      if (length(history_past$time) == 0) {
        cov_past[[cov_name]][i] <- 0
      } else {
        cov_past[[cov_name]][i] <- compute_covariates(
          CovariateName = cov_name,
          History       = history_past,
          Event         = e
        )
      }
      
      # current
      cov_current[[cov_name]][i] <- compute_covariates(
        CovariateName = cov_name,
        History       = history_current,
        Event         = e
      )
      
      # average
      cov_avg[[cov_name]][i] <- (cov_past[[cov_name]][i] + cov_current[[cov_name]][i]) / 2
    }
    
    # update history only when an event occurs 
    if (z_data[i] == 1) {
      history_data$time     <- c(history_data$time,     e$time)
      history_data$receiver <- c(history_data$receiver, e$receiver)
      history_data$sender   <- c(history_data$sender,   e$sender)
    }
  }
  
  # return the requested version
  if (type == "avg")     return(as.data.frame(cov_avg))
  if (type == "curr")    return(as.data.frame(cov_current))
  if (type == "past")    return(as.data.frame(cov_past))
}