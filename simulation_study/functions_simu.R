################################################################################
############################# HELPER FUNCTIONS #################################
################################################################################

# cross-reciprocity: exponentially decayed past events from the opposite network
# (reversed roles: r -> s)
compute_cross <- function(history_other, t, delta_rec){
  sapply(1:n, function(j){
    s <- dyads$sender[j]
    r <- dyads$receiver[j]
    times <- history_other$time[
      history_other$sender == r &
        history_other$receiver == s &
        history_other$time < t]
    sum(exp(-delta_rec * (t - times)))
  })
}

# within-repetition: exponentially decayed past events within the same network 
# (same direction: s -> r)
compute_rep_decay <- function(history, t, delta_rep){
  sapply(1:n, function(j){
    s <- dyads$sender[j]
    r <- dyads$receiver[j]
    times <- history$time[
      history$sender == s &
        history$receiver == r &
        history$time < t]
    sum(exp(-delta_rep * (t - times)))
  })
}

# construct case-control paired data: event dyad against random non-events for
# logistic regression
make_pairs <- function(events){
  lapply(events, function(ev){
    non_idx <- sample(setdiff(1:n, ev$dyad_idx), n_non_events)
    data.frame(
      time = ev$time,
      delta_cross = ev$cross_event - ev$cross_all[non_idx],
      delta_rep   = ev$rep_event   - ev$rep_all[non_idx],
      y = 1
    )
  }) |> dplyr::bind_rows()
}

# cross-repetition: covariate as exponentially decayed past events from the 
# opposite network (same direction: s -> r)
compute_cross_rep <- function(history_other, t, delta_rep){
  sapply(1:n, function(j){
    s <- dyads$sender[j]
    r <- dyads$receiver[j]
    times <- history_other$time[
      history_other$sender == s &
        history_other$receiver == r &
        history_other$time < t 
    ]
    sum(exp(-delta_rep * (t - times)))
  })
}

# construct case-control paired data for models with a single cross-effect covariate
make_pairs_cross <- function(events){
  if(length(events) == 0) return(NULL)
  lapply(events, function(ev){
    non_idx <- sample(setdiff(1:n, ev$dyad_idx), n_non_events)
    data.frame(
      time        = ev$time,
      delta_cross = ev$cross_event - ev$cross_all[non_idx],
      y           = 1
    )
  }) |> dplyr::bind_rows()
}
