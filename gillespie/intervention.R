library(ode2scm)

sde_sim <- function(transition_function, initial_states, times, interventions = NULL) {
  if(!is.null(interventions)) {
    for(int in names(interventions)){
      initial_states[[int]] <- interventions[[int]]
    }
  }
  print(initial_states)
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  t_delta <- times[2] - times[1]
  out <- as_tibble(
    smfsb::simTs(initial_states, times[1], times[length(times)], t_delta, transition_function)
  )
  out$time <- times
  out <- out[, c('time', setdiff(names(out), 'time'))]
  return(out)
}

gil <- function (N) 
{
  S = t(N$Post - N$Pre)
  v = ncol(S)

  return(function(x0, t0, deltat, ...) {
    t = t0
    x = x0
    print(x) # x is the initial state
    termt = t0 + deltat
    repeat {
      h = N$h(x, t, ...)
      print('h is...')
      print(h)
      h0 = sum(h)
      print('h0 is...')
      print(h0)
      if (h0 < 1e-10) t = 1e+99 else if (h0 > 1e+06) {
        t = 1e+99
        warning("Hazard too big - terminating simulation!")
      } else t = t + rexp(1, h0)
      print('t is..')
      print(t)
      if (t >= termt) return(x)
      j = sample(v, 1, prob = h)
      print('j is...')
      print(j)
      x = x + S[, j]
      print('let us see')
      print(x)
      print('------')

    }
    print('done repeating')
    print('------------------')
  })
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
  
  innerIntervention <- interventions
  
  sde$h <- function(states, t, parameters=rates, interventions = innerIntervention){
   #  print(interventions)
    # update the initial states
    if(!is.null(interventions)) {
      for(int in names(interventions)){
        states[[int]] <- interventions[[int]]
      }
    }
    print(states)
    with(as.list(c(states, parameters, interventions)), {
      # if intervene on PRaf, then both Raf and PRaf are fixed?
      # intervene on PPMek, then PMek and PPMek are fixed?
      
      #dE1 <- function(...) {
       # 0
      #}
      
      #dPRaf <- function(raf_activate, PRaf, E1, raf_deactivate, TRaf = 100, ...) {
       # raf_activate * (TRaf-PRaf) * E1 - raf_deactivate * PRaf
      #}
      
      #dPPMek <- function(mek_activate, mek_deactivate, PRaf, PPMek, TMek = 100, ...) {
       # (mek_activate ^ 2) * (PRaf ^ 2) * (TMek - PPMek) / mek_deactivate -       mek_activate * PRaf * PPMek -
        #  mek_deactivate * PPMek
      #}
      
      #dPPErk <- function(erk_activate, erk_deactivate, PPMek, PPErk, TErk = 100, ...) {
       # (erk_activate ^ 2) * (PPMek ^ 2) * (TErk - PPErk) / erk_deactivate -      erk_activate * PPMek * PPErk -
        #  erk_deactivate * PPErk
      #}
      
      out <- c(
        raf_activate * Raf * E1,
        raf_deactivate * PRaf,
        0,
        0,
        0,
        0,
        erk_activate * PPMek * Erk,
        erk_deactivate * PErk,
        erk_activate * PPMek * PErk,
        erk_deactivate * PPErk
      )
      #out[1] <- 0
      #out[2] <- 0   
      if (!is.null(interventions$Raf) || !is.null(interventions$PRaf)) {
        print('intervene on raf')
        out[1] <- 0
        out[2] <- 0        
      } else if (!is.null(interventions$Mek) || 
                 !is.null(interventions$PMek) || 
                 !is.null(interventions$PPMek)) {
        print('intervene on mek')
        out[3] <- 0
        out[4] <- 0
        out[5] <- 0
        out[6] <- 0
      } 
      #else {
        # print('intervene on erk')
        #out[7] <- 0
        #out[8] <- 0
        #out[9] <- 0
        #out[10] <- 0
      #}

      return(out)
    })
  }
  transition_function <- StepGillespie(sde)
  # h = N$h(x, t, ...)
  # h0 = sum(h)
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

initial_states <-  list(E1=1, Raf=100, PRaf=0, Mek=50, PMek=20, PPMek=30, Erk=100, PErk=0, PPErk=0)

intervention_raf <- list(Raf = 70, PRaf = 30)
intervention_mek <- list(Mek=50, PMek=20, PPMek=30)

times <- seq(0, 0.5, by = .1)
faster_rates <- lapply(rates, `*`, 20) # WTF
stoc_transition_func <- mapk_sde(initial_states, faster_rates
                                 )
sde_out <- sde_sim(stoc_transition_func, initial_states, times
                   )

sde_out[nrow(sde_out),]$Raf
sde_out[nrow(sde_out),]$PRaf
sde_out[nrow(sde_out),]$Mek
sde_out[nrow(sde_out),]$PMek
sde_out[nrow(sde_out),]$PPMek
sde_out[nrow(sde_out),]$Erk
sde_out[nrow(sde_out),]$PErk
sde_out[nrow(sde_out),]$PPErk




