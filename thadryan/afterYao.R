suppressWarnings(suppressMessages(library(gdata)))
suppressWarnings(suppressMessages(library(ode2scm)))

ode_sim <- function(transition_function, initial_states, times){
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  # print(initial_states)
  # print(transition_function)
  rates <- attr(transition_function, 'rates')
  rates <- structure(as.numeric(rates), names = names(rates))
  # print(rates)
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
  innerinterventions <- interventions #(without mapkDo)
  transition_function <- function(t,
                                  states = innerStates,
                                  rates = innerRates,
                                  interventions = innerinterventions) {
    for(int in names(interventions)){
      states[[int]] <- interventions[[int]]
    }

    states <- as.list(states)
    rates <- as.list(rates)

    # print('check variables')
    # print('PRaf')
    # print(states$PRaf)
    # print('PPMek')
    # print(states$PPMek)
    # print('PPErk')
    # print(states$PPErk)


    dE1 <- 0
    if(!is.null(interventions$PRaf)) {
      dPErk <- 0
    } else{
      dPRaf <- rates$raf_activate * (100-states$PRaf) * states$E1 -
        rates$raf_deactivate * states$PRaf
    }
    if(!is.null(interventions$PPMek)){
      dPPMek <- 0
    } else{
      dPPMek <- (rates$mek_activate ^ 2) * (states$PRaf ^ 2) * (100 - states$PPMek) /
        rates$mek_deactivate - rates$mek_activate * states$PRaf * states$PPMek -
        rates$mek_deactivate * states$PPMek
    }
    if(!is.null(interventions$PPErk)){
      dPPErk <- 0
    }else{
      dPPErk <- (rates$erk_activate ^ 2) * (states$PPMek ^ 2) * (100 - states$PPErk) /
        rates$erk_deactivate - rates$erk_activate * states$PPMek * states$PPErk -
        rates$erk_deactivate * states$PPErk
    }

    # print('dPRaf')
    # print(dPRaf)
    # print('dPPMek')
    # print(dPPMek)
    # print('dPPErk')
    # print(dPPErk)

    list(c(dE1, dPRaf, dPPMek, dPPErk))
  }
  attr(transition_function, 'rates') <- as.list(rates)
  attr(transition_function, 'states') <- as.list(states)
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

# so I guess we don't need Raf, Mek, PMek, Erk and PErk values!
initial_states <-  list(
  E1 = 1,
  PRaf = 0,
  PPMek = 20,
  PPErk = 0
)

times <- seq(0, 120, by = .1)

# mapk_ode_instance <- mapk_ode(initial_states, rates)
# mapk_ode_intervention <- do(mapk_ode_instance, list(PPMek = k))
# ode_out <- ode_sim(mapk_ode_intervention, initial_states, times)

## with intervention
mapkDo <- function(model, doIntervention)
{
  model(interventions = doIntervention)
}

## without intervention
print("Without Intervention")
firstModel <- mapk_ode_robert(initial_states, rates)
firstModel()
#testOut <- ode_sim(firstModel, initial_states, times)

# # result
# PRafM <- testOut[nrow(testOut),]$PRaf
# PPMekM <- testOut[nrow(testOut),]$PPMek
# PPErkM <- testOut[nrow(testOut),]$PPErk
#
# print("With Intervention not working")
# intervened <- mapkDo(firstModel, list(PPMek=10))
# intervened()
# testOutIntervention <- ode_sim(intervened, initial_states, times)
#
# print("With Intervention working")
# intervented2 <- mapk_ode_robert(initial_states, rates, list(PPMek=20)) # have to change initial_states
# intervented2()
# testOutIntervention2 <- ode_sim(intervented2, initial_states, times)
#
# PRafM_intervened <- testOutIntervention2[nrow(testOutIntervention2),]$PRaf
# PPMekM_intervened <- testOutIntervention2[nrow(testOutIntervention2),]$PPMek
# PPErkM_intervened <- testOutIntervention2[nrow(testOutIntervention2),]$PPErk
