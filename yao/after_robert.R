library(gdata)
library(ode2scm)
ode_sim <- function(transition_function, initial_states, times){
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  rates <- attr(transition_function, 'rates')
  rates <- structure(as.numeric(rates), names = names(rates))
  as_tibble(
    as.data.frame(
      deSolve::ode(
        y = initial_states,
        times = times,
        func = transition_function,
        parms = rates
      )
    )
  )
}

mapk_ode_simplified <- function(states, rates) {
  transition_function <- function(t, interventions=NULL) {
    if(!is.null(interventions)) {
      states = update.list(states, interventions)
    }
    with(as.list(c(innerStates, innerRates)), {
      dE1 <- 0
      dPRaf <- raf_activate * Raf * E1 - raf_deactivate * PRaf
      dPPMek <- (mek_activate ^ 2) * (PRaf ^ 2) * (Mek + PMek) /
        mek_deactivate - mek_activate * PRaf * PPMek - mek_deactivate * PPMek
      dPPErk <- (erk_activate ^ 2) * (PPMek ^ 2) * (Erk + PErk) /
        erk_deactivate - erk_activate * PPMek * PPErk - erk_deactivate * PPErk
    })
  }
}
if(!is.null(interventions)) {
  states = update.list(states, interventions)
}



mapk_ode_robert <- function(states, rates, interventions=NULL) {
  innerRates <- rates
  innerStates <- states
  transition_function <- function(t, states = innerStates, rates = innerRates,
                                  interventions = NULL) {
    for(int in names(interventions)){
      states[[int]] <- interventions[[int]]
    }
    
    states <- as.list(states)
    rates <- as.list(rates)
    
    dE1 <- 0
    if(!is.null(interventions$PRaf)) {
      dPErk <- 0
    } else{
      dPRaf <- rates$raf_activate * (100-states$PRaf) * states$E1 - rates$raf_deactivate * states$PRaf
    }
    if(!is.null(interventions$PPMek)){
      dPPMek <- 0
    } else{
      dPPMek <- (rates$mek_activate ^ 2) * (states$PRaf ^ 2) * (100 - states$PPMek) /
        rates$mek_deactivate - rates$mek_activate * states$PRaf * states$PPMek - rates$mek_deactivate * states$PPMek
    }
    if(!is.null(interventions$PPErk)){
      dPPErk <- 0
    }else{
      dPPErk <- (rates$erk_activate ^ 2) * (states$PPMek ^ 2) * (100 - states$PPErk) /
        rates$erk_deactivate - rates$erk_activate * states$PPMek * states$PPErk - rates$erk_deactivate * states$PPErk
    }
    
    list(c(dE1, dPRaf, dPPMek, dPPErk))
  }
  attr(transition_function, 'rates') <- as.list(innerRates)
  attr(transition_function, 'states') <- as.list(innerStates)
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
  PRaf = 20,
  PPMek = 0,
  PPErk = 0
)

times <- seq(0, 120, by = .1)

# mapk_ode_instance <- mapk_ode(initial_states, rates)
# mapk_ode_intervention <- do(mapk_ode_instance, list(PPMek = k))
# ode_out <- ode_sim(mapk_ode_intervention, initial_states, times)

testModelInstance <- mapk_ode_robert(initial_states, rates)
testOut <- ode_sim(testModelInstance, initial_states, times)

# result
PRafM <- testOut[nrow(testOut),]$PRaf
PPMekM <- testOut[nrow(testOut),]$PPMek
PPErkM <- testOut[nrow(testOut),]$PPErk

# so I guess we don't need Raf, Mek, PMek, Erk and PErk values!
initial_states <-  list(
  E1 = 1,
  PRaf = 20,
  PPMek = 0,
  PPErk = 0
)

mapkDo <- function(model, doIntervention)
{
  intervened_model = function() {
    model(interventions = doIntervention)
  }
  return (intervened_model)
}

print("Before...")
firstModel <- mapk_ode_robert(initial_states, rates)
firstModel()

print("After...")
intervened <- mapkDo(firstModel, list(PPMek=10))
intervened()
