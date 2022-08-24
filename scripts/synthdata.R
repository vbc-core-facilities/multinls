#!/usr/bin/env Rscript

# This script generates synthetic data
# for the FRAP data analysis example.
# András Aszódi, 2022-08-24

# The nonlinear function F(t) 
source("frapfun.R")

# This function generates synthetic data for FRAP nonlinear fit analysis.
# par.df = Data frame of the "known" parameters:
#   the columns are thalf, f0, finf
#   the rows correspond to the parameters in one group.
#   There will be as many observation groups as rows in this data frame.
# tmax.fact = The range of time points from 0 up to this value times the
#   maximum thalf in `par.df`. Default is 8, no need to change.
# nobs = The number of observations (time points)
# noise.sd = If not 0.0, then a Normal random number with this standard deviation
#   will be used to add a relative error to the observations.
# Returns a data frame, the first column is called "times" 
#   and contains the time points (the independent variables)
#   the rest of the columns contain the measurements at the time points
#   and are named "A", "B", "C", ... etc.
synth.data <- function(par.df, tmax.fact=8.0, nobs=50, noise.sd=0.0) {
  # time runs from 0.0 to 8 times the maximal `thalf`
  times <- seq(0.0, tmax.fact*max(par.df$thalf), length.out=nobs)
  
  # save everything with this many significant digits
  DIGITS <- 4
  
  # construct observations for each row in `par.df` (different curves)
  # a matrix will be returned
  dm <- apply(par.df, MARGIN=1, function(pars) {
    # `pars` is a vector, give it back the `par.df` column names and convert to list
    names(pars) <- colnames(par.df)
    plist <- as.list(pars)
    # `nobs`-length vector of noiseless response variable
    # for a fixed set of parameters at all time points
    fts <- sapply(times, function(tm) {
      ft <- do.call(frap.fun, c(t=tm, plist)) # concatenate to form parameter list for invocation
      return(ft)
    })
    # add some relative noise
    noise <- rnorm(nobs, sd=noise.sd)
    return(signif(fts*(1.0 + noise), DIGITS))
  })
  
  # prepend the times as the first column
  dm <- cbind(signif(times, DIGITS), dm)
  
  # Convert the matrix to a data frame
  # name the response variable columns with capital letters
  # hope that ncol <= 26 :-)
  colnames(dm) <- c("times", LETTERS[1:(ncol(dm)-1)])
  d <- as.data.frame(dm)
  return(d)
}

# == MAIN ==

# The "known" parameters
# there are 3 parameter sets corresponding to the rows of this data frame
par.df <- data.frame(thalf=c(11.0, 12.0, 14.0), f0=c(0.1, 0.3, 0.2), finf=c(2.1, 3.7, 5.8))
row.names(par.df) <- c("A","B","C")
cat("Table of the 'known' parameters for the three groups:\n")
print(par.df)

# Synthetic dataset with a little noise
set.seed(137) # reproducible
d <- synth.data(par.df, noise.sd=0.01)
# save
write.csv(d, file="../data/synthdata.csv", row.names = F)
