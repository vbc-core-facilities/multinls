# Source, not execute.

# This script contains the R implementation
# of the nonlinear function describing the photobleaching time course
# in a FRAP experiment.
# Based on Equation 12 (page 70) from the paper
# Yguerabide J. et al, Biophys. J. 39:69-75 (1982).

# t is time, the independent variable.
# thalf, f0, finf are the three parameters
# corresponding to t_{1/2}, F(0) and F(\infty) (in LaTeX notation) in the paper.
# Returns F(t)
frap.fun <- function(t, thalf, f0, finf) {
  trel <- t/thalf
  ft <- (f0 + finf*trel)/(1.0 + trel)
  return(ft)
}
