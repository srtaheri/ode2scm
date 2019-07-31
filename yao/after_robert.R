suppressWarnings(suppressMessages(library(gdata)))
suppressMessages(suppressMessages(library(tibble)))

ode_sim <- function(transition_function, initial_states, times){
  print("called")
  initial_states <- structure(as.numeric(initial_states), names = names(initial_states))
  rates <- attr(transition_function, 'rates')
  rates <- structure(as.numeric(rates), names = names(rates))
  rates <- as.list(rates)
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

mapk_ode_robert2 <- function(states, rates, interventions = NULL) {

  # avoid recusive argument error
  innerRates <- rates
  innerStates <- states

  transition_function <- function(t,
    states = innerStates, rates = innerRates) {

    if(!is.null(interventions)) {
      for(int in names(interventions)){
        states[[int]] <- interventions[[int]]
      }
    }

    states <- as.list(states)
    rates <- as.list(rates)

    dE1 <- 0

    if(!is.null(interventions$PRaf)) {
      dPRaf <- 0
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
    print("My states are:")
    print(unlist(states))
    list(c(dE1, dPRaf, dPPMek, dPPErk))
  }
  attr(transition_function, 'rates') <- as.list(rates)
  attr(transition_function, 'states') <- as.list(states)
  return(transition_function)
}
model = function(states, rates) {
  return(list("dE1" = 0,
              "dPRaf" = rates$raf_activate * (100-states$PRaf) * states$E1 -
                rates$raf_deactivate * states$PRaf,
              "dPPMek" = (rates$mek_activate ^ 2) * (states$PRaf ^ 2) * (100 - states$PPMek) /
                rates$mek_deactivate - rates$mek_activate * states$PRaf * states$PPMek -
                rates$mek_deactivate * states$PPMek,
              "dPPErk" = (rates$erk_activate ^ 2) * (states$PPMek ^ 2) * (100 - states$PPErk) /
                rates$erk_deactivate - rates$erk_activate * states$PPMek * states$PPErk -
                rates$erk_deactivate * states$PPErk)
  )
}

mapk_ode_robert3 <- function(model, states, rates, interventions = NULL) {
  
  # avoid recusive argument error
  innerRates <- rates
  innerStates <- states
  
  transition_function <- function(t,
                                  states = innerStates,
                                  rates = innerRates) {
    
    if(!is.null(interventions)) {
      for(int in names(interventions)){
        states[[int]] <- interventions[[int]]
      }
    }
    
    states <- as.list(states)
    rates <- as.list(rates)
    
    d_intervene <- list()

    
    dE1 <- 0
    
    dPRaf <- rates$raf_activate * (100-states$PRaf) * states$E1 -
      rates$raf_deactivate * states$PRaf
  
    dPPMek <- (rates$mek_activate ^ 2) * (states$PRaf ^ 2) * (100 - states$PPMek) /
      rates$mek_deactivate - rates$mek_activate * states$PRaf * states$PPMek -
      rates$mek_deactivate * states$PPMek
    
    dPPErk <- (rates$erk_activate ^ 2) * (states$PPMek ^ 2) * (100 - states$PPErk) /
      rates$erk_deactivate - rates$erk_activate * states$PPMek * states$PPErk -
      rates$erk_deactivate * states$PPErk
    
    for (elem in names(interventions)) {
      d_intervene <- paste0("d", elem, " <- 0")
      eval(parse(text=d_intervene))
    }
    
    print("My states are:")
    print(unlist(states))
    list(c(dE1, dPRaf, dPPMek, dPPErk))
    print("dPRaf")
    print(dPRaf)
    print("dPPMek")
    print(dPPMek)
    print("dPPErk")
    print(dPPErk)
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
initial_states <- list(
  E1 = 1,
  PRaf = 0,
  PPMek = 0,
  PPErk = 0
)

times <- seq(0, 120, by = .1)

firstModel <- mapk_ode_robert2(initial_states, rates)
firstModel()
ode_out <- ode_sim(firstModel, initial_states, times)

# I am not sure there is really a point to making do a function
# if we're using a deterministic functtion
secondModel <- mapk_ode_robert2(initial_states, rates,
                                interventions = list(PPMek = 30))
secondModel()

ode_out_intervention <- ode_sim(secondModel, initial_states, times)
# # result
# problem: intervened variable will show up as 0
PRafM <- ode_out_intervention[nrow(ode_out_intervention),]$PRaf
PPMekM <- ode_out_intervention[nrow(ode_out_intervention),]$PPMek
PPErkM <- ode_out_intervention[nrow(ode_out_intervention),]$PPErk

test_intervene <- list(PRaf = 10)
d_intervene <- list()
for (elem in names(test_intervene)) {
  d_intervene <- c(d_intervene, paste0("d", elem, " <- 0"))
}


thirdModel <- mapk_ode_robert3(initial_states, rates,
                                interventions = list(PPMek = 30))
thirdModel()
