library(gdata)

# need to modify for interventions: in case of intervention, do not calculate the equilibrium value? 
mapk_ode_equilM <- function(states, rates, interventions = NULL) {
  innerState <- states
  innerRate <- rates
  transition_f <- function(states = innerState, rates = innerRate) {

    with(as.list(c(states, rates)), {
      w3 <- raf_activate/raf_deactivate
      w2 <- mek_activate/mek_deactivate
      w1 <- erk_activate/erk_deactivate
      
      # can we assume that t1, t2 and t3 are all 100? b/c in initial_states we're only keeping PRaf, PPMek and PPErk
      t3 <- Raf + PRaf
      t2 <- Mek + PMek + PPMek
      t1 <- Erk + PErk + PPErk
      
      if(!is.null(interventions$PRaf)) {
        k3 <- interventions$PRaf
      } else {
        u3 <- w3 * E1
        k3 <- t3 * (u3/(1+u3))      
      }
      
      if(!is.null(interventions$PPMek)) {
        k2 <- interventions$PPMek
      } else {
        u2 <- w2 * k3
        k2 <- t2 * ((u2^2)/(1 + u2 + u2^2))
      }
  
      if(!is.null(interventions$PPErk)) {
        k1 <- interventions$PPErk
      } else {
        u1 <- w1 * k2
        k1 <- t3 * ((u1^2)/(1 + u1 + u1^2))      
      }
  
      list(k3,k2,k1)
      
      })
    }
  return(transition_f)
  
}
intervention <- list(PPErk = 30)

rates <- list(
  raf_activate = 0.1,
  raf_deactivate = 0.1,
  mek_activate = 0.1,
  mek_deactivate = 2.0,
  erk_activate = 0.1,
  erk_deactivate = 1.0
)

initial_states <-  list(
  E1 = 1,
  Raf = 100,
  PRaf = 0,
  Mek = 100,
  PMek = 0,
  PPMek = 0,
  Erk = 100,
  PErk = 0,
  PPErk = 0
)

mapk_ode_equilM(initial_states, rates, intervention)()