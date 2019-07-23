diff_equations = list(
  E1 =0,
  Raf =-raf_activate * Raf * E1 + raf_deactivate * PRaf,
  PRaf =raf_activate * Raf * E1 - raf_deactivate * PRaf,
  Mek =-mek_activate * Mek * PRaf + mek_deactivate * PMek,
  PMek =mek_activate * Mek * PRaf +
    mek_deactivate * PPMek -
    mek_deactivate * PMek -
    mek_activate * PMek * PRaf,
  PPMek =mek_activate * PMek * PRaf - mek_deactivate * PPMek,
  Erk =-erk_activate * Erk * PPMek + erk_deactivate * PErk,
  PErk =erk_activate * Erk * PPMek +
    erk_deactivate * PPErk -
    erk_deactivate * PErk -
    erk_activate * PErk * PPMek,
  PPErk =erk_activate * PErk * PPMek - erk_deactivate * PPErk
)

x = NULL
x = paste0(x, "{")
x = paste0(x, "print('hello world')")
x = paste0(x, "}")

eval(parse(text=x))


{
  dE1 <- 0
  dRaf <- -raf_activate * Raf * E1 + raf_deactivate * PRaf
  dPRaf <- raf_activate * Raf * E1 - raf_deactivate * PRaf
  dMek <- -mek_activate * Mek * PRaf + mek_deactivate * PMek
  dPMek <- mek_activate * Mek * PRaf +
    mek_deactivate * PPMek -
    mek_deactivate * PMek -
    mek_activate * PMek * PRaf
  dPPMek <- mek_activate * PMek * PRaf - mek_deactivate * PPMek
  dErk <- -erk_activate * Erk * PPMek + erk_deactivate * PErk
  dPErk <- erk_activate * Erk * PPMek +
    erk_deactivate * PPErk -
    erk_deactivate * PErk -
    erk_activate * PErk * PPMek
  dPPErk <- erk_activate * PErk * PPMek - erk_deactivate * PPErk
  
  list(c(dE1, dRaf, dPRaf, dMek, dPMek, dPPMek, dErk, dPErk, dPPErk))
}


# active Raf = PRAF
# active Mek = PPMek
# active Erk = PPErk

mapk_ode <- function(states, rates, interventions){
  transition_function <- function(t, states, rates) {
    with(as.list(c(states, rates)), {
      dE1 <- 0
      dRaf <- -raf_activate * Raf * E1 + raf_deactivate * PRaf
      dPRaf <- raf_activate * Raf * E1 - raf_deactivate * PRaf
      dMek <- -mek_activate * Mek * PRaf + mek_deactivate * PMek
      dPMek <- mek_activate * Mek * PRaf +
        mek_deactivate * PPMek -
        mek_deactivate * PMek -
        mek_activate * PMek * PRaf
      dPPMek <- mek_activate * PMek * PRaf - mek_deactivate * PPMek
      dErk <- -erk_activate * Erk * PPMek + erk_deactivate * PErk
      dPErk <- erk_activate * Erk * PPMek +
        erk_deactivate * PPErk -
        erk_deactivate * PErk -
        erk_activate * PErk * PPMek
      dPPErk <- erk_activate * PErk * PPMek - erk_deactivate * PPErk
      
      list(c(dE1, dRaf, dPRaf, dMek, dPMek, dPPMek, dErk, dPErk, dPPErk))
    })
  }
  attr(transition_function, 'rates') <- rates
  return(transition_function)
}