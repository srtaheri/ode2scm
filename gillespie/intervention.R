library(ode2scm)

sde_sim <- function(transition_function, initial_states, times, interventions = NULL) {
  if(!is.null(interventions)) {
    for(int in names(interventions)){
      initial_states[[int]] <- interventions[[int]]
    }
  }
  #print(initial_states)
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  t_delta <- times[2] - times[1]
  out <- as_tibble(
    smfsb::simTs(initial_states, times[1], times[length(times)], t_delta, transition_function)
  )
  out$time <- times
  out <- out[, c('time', setdiff(names(out), 'time'))]
  return(out)
}


mapk_sde <- function(states, rates, interventions = NULL){
  sde <- list()
  
  sde$Pre <- matrix(
    c(
      1, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1
    ), nrow=10, ncol=9, byrow=T
  )
  colnames(sde$Pre) <- c("E1","Raf", "PRaf", "Mek", "PMek", "PPMek", "Erk", "PErk", "PPErk")
  sde$Post <- matrix(
    c(
      1, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0
    ), nrow=10, ncol=9, byrow=T
  )
  colnames(sde$Post) <- c("E1","Raf", "PRaf", "Mek", "PMek", "PPMek", "Erk", "PErk", "PPErk")
  innerIntervention <- interventions
  
  sde$h <- function(states, t, parameters=rates, interventions = innerIntervention){
   #  print(interventions)
    # update the initial states
    if(!is.null(interventions)) {
      for(int in names(interventions)){
        states[[int]] <- interventions[[int]]
      }
    }
    #print(states)
    with(as.list(c(states, parameters, interventions)), {
      if(!is.null(interventions)) {
        for(int in names(interventions)){
          sde$Pre[,int] <- 0
          sde$Post[,int] <- 0
        }
      }
      out <- c(
        raf_activate * Raf * E1,
        raf_deactivate * PRaf,
        mek_activate * PRaf * Mek,
        mek_deactivate * PMek,
        mek_activate * PRaf * PMek,
        mek_deactivate * PPMek,
        erk_activate * PPMek * Erk,
        erk_deactivate * PErk,
        erk_activate * PPMek * PErk,
        erk_deactivate * PPErk
      )
      print(out)

      return(out)
    })
  }
  transition_function <- StepGillespie(sde)
  return(transition_function)
}

rates <- list(
  raf_activate=0.1,
  raf_deactivate=0.1,
  mek_activate=0.1,
  mek_deactivate=2.0,
  erk_activate=0.1,
  erk_deactivate=1.0
)

initial_states <-  list(E1=1, Raf=100, PRaf=0, Mek=100, PMek=0, PPMek=0, Erk=100, PErk=0, PPErk=0)

intervention_raf <- list(Raf = 70, PRaf = 30)
intervention_mek <- list(Mek=50, PMek=20, PPMek=30)
intervention_mek_raf <- list(Raf = 70, PRaf = 30,Mek=50, PMek=20, PPMek=30)

times <- seq(0, 0.5, by = .1)
faster_rates <- lapply(rates, `*`, 20) # WTF
stoc_transition_func <- mapk_sde(initial_states, faster_rates,intervention_mek_raf
                                 )
sde_out <- sde_sim(stoc_transition_func, initial_states, times,intervention_mek_raf
                   )

sde_out[nrow(sde_out),]$Raf
sde_out[nrow(sde_out),]$PRaf
sde_out[nrow(sde_out),]$Mek
sde_out[nrow(sde_out),]$PMek
sde_out[nrow(sde_out),]$PPMek
sde_out[nrow(sde_out),]$Erk
sde_out[nrow(sde_out),]$PErk
sde_out[nrow(sde_out),]$PPErk




