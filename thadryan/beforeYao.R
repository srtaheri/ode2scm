suppressWarnings(suppressMessages(library(gdata)))
suppressMessages(suppressMessages(library(tibble)))

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

mapk_ode_robert <- function(states, rates){

  # I don't like this any more than you do \_(> _ <)_/
  # ...but we get a recusive argument error if we don't do this...
  innerRates <- rates
  innerStates <- states

  # ...because these need defaults that have different names
  transition_function <- function(t, states = innerStates, rates = innerRates,
                                  interventions = NULL) {

    if(!is.null(interventions)) {
      print(paste("interventing on: ", names(interventions),
        "value of", unlist(interventions)))
    }

    for(int in names(interventions)){
      states[[int]] <- interventions[[int]]
    }

    dE1 <- 0
    if(!is.null(interventions$PRaf)) {
      dPErk <- 0
    } else {
      dPRaf <- rates$raf_activate * (100-states$PRaf) *
        states$E1 - rates$raf_deactivate * states$PRaf
    }


    if(!is.null(interventions$PPMek)) {
      dPPMek <- 0
    } else {
      dPPMek <- (rates$mek_activate ^ 2) * (states$PRaf ^ 2) * 100 /
        rates$mek_deactivate - rates$mek_activate * states$PRaf *
          states$PPMek - rates$mek_deactivate * states$PPMek
    }


    if(!is.null(interventions$PPErk)){
      dPPErk <- 0
    } else {
      dPPErk <- (rates$erk_activate ^ 2) * (states$PPMek ^ 2) * 100 /
        rates$erk_deactivate - rates$erk_activate * states$PPMek *
        states$PPErk - rates$erk_deactivate * states$PPErk
    }

    list(dE1, dPRaf, dPPMek, dPPErk)
  }

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
  Raf = 80,
  PRaf = 20,
  Mek = 100,
  PMek = 0,
  PPMek = 0,
  Erk = 100,
  PErk = 0,
  PPErk = 0
)

times <- seq(0, 120, by = .1)

mapkDo <- function(model, doIntervention)
{
  # if we intervene here we're going to get a new function
  x <- model(interventions = doIntervention)
  list(fcn = model, data = x)
}

print("Before...")
firstModel <- mapk_ode_robert(initial_states, rates)
firstModel()

print("After...")
intervened <- mapkDo(firstModel, list(PPMek = 10))
# why no ()?
intervened

f <- getElement(intervened, "fcn")
print(class(f))
ode_out <- ode_sim(f, initial_states, times)
print("what is your problem")
