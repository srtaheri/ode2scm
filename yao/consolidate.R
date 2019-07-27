library(gdata)

mapk_ode_equilM <- function(states, rates) {
  transition_f <- function() {
    w3 <- raf_activate/raf_deactivate
    w2 <- mek_activate/mek_deactivate
    w1 <- erk_activate/erk_deactivate
    t3 <- Raf + PRaf
    t2 <- Mek + PMek + PPMek
    t1 <- Erk + PErk + PPErk
    
    u3 <- w3 * E1
    k3 <- t3 * (u3/(1+u3))
    
    u2 <- w2 * k3
    k2 <- t2 * ((u2^2)/(1 + u2 + u2^2))
    
    u1 <- w1 * k2
    k1 <- t3 * ((u1^2)/(1 + u1 + u1^2))
    
    list(k3,k2,k1)
    print(list(k3,k2,k1))
  }
  return(transition_f)
}


mapk_ode_simplified <- function(states, rates) {
  
  transition_function <- function(t, interventions=NULL) {
    
    # intervene if we're supposed to
    if(!is.null(interventions)) {
      states = update.list(states, interventions)
    }
    
    # run the code normally
    with(as.list(c(states, rates)),{
      
      print("My states are...")
      print(unlist(states))
      
      dE1 <- 0
      #print("PPMek intervention is:")
      #print(is.null(interventions$PRaf))
      if(is.null(interventions$PRaf)) {
        dPRaf <- raf_activate * Raf * E1 - raf_deactivate * PRaf
      }else {
        dPRaf = 0
      }
      if(is.null(interventions$PPMek)) {
        dPPMek <- (mek_activate ^ 2) * (PRaf ^ 2) * (Mek + PMek) /
          mek_deactivate - mek_activate * PRaf * PPMek - mek_deactivate * PPMek
      }else {
        dPPMek = 0
      }
      if(is.null(interventions$PPErk)) {
        dPPErk <- (erk_activate ^ 2) * (PPMek ^ 2) * (Erk + PErk) /
          erk_deactivate - erk_activate * PPMek * PPErk - erk_deactivate * PPErk
      }else {
        dPPErk = 0
      }
      list(c(dE1, dPRaf, dPPMek, dPPErk))
    })
  }
  attr(transition_function, 'rates') <- rates
  return(transition_function)
}

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
  Raf = 80,
  PRaf = 20,
  Mek = 100,
  PMek = 0,
  PPMek = 0,
  Erk = 100,
  PErk = 0,
  PPErk = 0
)

mapkDo <- function(model, intervention)
{
  model(interventions = intervention)
}

print("Before...")
testOut <- mapk_ode_simplified(initial_states, rates)
testOut()

print("After...")
new <- mapkDo(testOut, list(PPMek = 10))
