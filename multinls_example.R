#!/usr/bin/env Rscript

# Example script for FRAP analysis
# Runs a (multiple) nonlinear regression.
# For Stela Jelenic, Saha/IMBA.
# 2022-08-23 András Aszódi

# Literature:
# Yguerabide J. et al, Biophys. J. 39:69-75 (1982)

# This package provides the `melt` function
library("reshape2")

source("scripts/frapfun.R")
source("scripts/nlsplotte.R")

# -- Functions --

# Get initial rough estimates for the parameters.
# d = the original data frame with columns "times", "A", "B", "C", ...
init.params <- function(d) {
  # calculate the average of all groups for each time point
  # (the first column, "times", is omitted)
  fmeans <- apply(d[, -1], MARGIN=1, mean)
  
  # the first value
  f0 <- fmeans[1]
  
  # the last value
  nobs <- length(fmeans)  # number of observations
  finf <- fmeans[nobs]
  
  # thalf corresponds to the time where F is approximately "halfway" between F0 and Finf
  # we _assume_ that the average f(t) values increment nicely, this is not foolproof
  fhalf <- (finf + f0)/2
  thalf <- (d$times[nobs] + d$times[1])/2 # simplest estimate: halfway between tmin and tmax
  tprev <- d$times[1]
  fprev <- fmeans[1]
  for(r in 2:nobs) {
    tcur <- d$times[r]
    fcur <- fmeans[r]
    if (fprev <= fhalf && fhalf <= fcur) {
      thalf <- (tcur + tprev)/2.0
      break
    }
    tprev <- tcur
    fprev <- fcur
  }
  
  # return the estimates in a list
  return(list(thalf=thalf, f0=f0, finf=finf))
}

# Convert the "group-wise fit" coefficient vector to a data frame
# with a similar layout as the known parameters `par.df`
convert.coeffs <- function(coefs) {
  c.mat <- matrix(coefs, ncol=3)
  # the column order is "thalf", "f0", "finf"
  c.df <- as.data.frame(c.mat)
  colnames(c.df) <- c("thalf", "f0", "finf")
  rownames(c.df) <- NULL
  return(c.df)
}

# -- Convenience functions --

# Shows the head and the tail of a data frame.
head.tail <- function(df) {
  cat("First few lines:\n")
  print(head(df))
  cat("Last few lines:\n")
  print(tail(df))
}

# == MAIN ==

# Read the synthetic data set.
# Look at "scripts/synthdata.R" for the details how they were generated.
d <- read.csv("data/synthdata.csv")
cat("Original data set:\n")
head.tail(d)

# "tidy" the data.
# In the original data frame `d` each group has a column "A", "B", "C", ...
# The fits and the plots require that these are rearranged
# into a new data frame with the following columns:
# "times" are the abscissae, repeated for each group
# "grp" contains the labels "A","B","C" for the groups
# "ft" contains the corresponding F(t) values
ds <- melt(d, id.vars="times", variable.name="grp", value.name="ft")
cat("Rearranged (tidied) data set\n")
head.tail(ds)

# Starting parameters: rough estimates
start.pars <- init.params(d)
cat("Initial parameters, rough estimate:\n")
print(start.pars)

# NLS fit to the pooled observations
# This is our "baseline". If we get a "good fit",
# then the groups are not very different.
pooled.fit <- nls(ft ~ frap.fun(times, thalf, f0, finf), data=ds,
                  start=start.pars)
print(summary(pooled.fit))

# NLS fit to each of the groups
# the starting parameters must be "replicated" as many times as there are groups
ngroup <- ncol(d)-1
multi.start.pars <- lapply(start.pars, function(p) { rep(p, ngroup)})

# Fit the function to all groups separately.
# `nls` does this by "indexing" the parameters: in the formula (first parameter)
# note the `thalf[grp]` etc. expressions
# this is not explained in the help for `nls`, only shown as the very last example
group.fit <- nls(ft ~ frap.fun(times, thalf[grp], f0[grp], finf[grp]), data=ds,
               start=multi.start.pars)
print(summary(group.fit))

# Compare the two fits using an "ANOVA" F-test
# If this test result is "significant", it means that the group-wise fits
# describe the data "better" than the pooled fit, i.e. the group-wise
# estimated parameters differ.
print(anova(pooled.fit, group.fit))

# Organise the group-wise fit coefficients into a data frame
# the columns are "thalf", "f0", "finf"
# the rows correspond to the fitted parameters in each group.
#
# FYI the "true" parameters were (see "scripts/synthdata.R"):
#   thalf  f0 finf
# A    11 0.1  2.1
# B    12 0.3  3.7
# C    14 0.2  5.8
#
c.df <- convert.coeffs(coef(group.fit))
cat("Estimated parameter table:\n")
print(c.df)

## ggplot the observations and the 3 fitted functions
pf <- nls.plot(ds, c.df)
save.png("plots/fits", pf)
cat("Look at the plot 'plots/fits.png' :-)\n")

