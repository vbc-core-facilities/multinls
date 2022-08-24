# Source, not execute

# Plotting functions.

library("ggplot2")

# Saves a GG plot in PNG format.
# Make sure the script that `source`-s this also loads the "ggplot2" package.
# basename = the name of the file without the ".png" extension which will be added automatically
save.png <- function(basename, plot=last_plot(), w=12, h=8) {
  filename <- paste(basename, "png", sep=".")
  ggsave(filename, plot=plot, width=w, height=h)
}

# Ad-hoc data plotter.
# This function plots the observations as dots and optionally
# the fitted functions as black lines over the groups.
# ds: a data frame with three columns: "times","grp","ft"
#     where "times" is the independent (predictor) variable, the time points;
#     "grp" is a factor indicating the groups that are fitted separately, like "A","B","C", ...
#     "ft" is the dependent (response) variable F(t), the observed fluorescence values.
# If `c.df` is not NULL, then it must be a data frame
#     containing the fitted parameters. The columns are "thalf", "f0", "finf"
#     corresponding to the parameters of the function `frap.fun` (see "frapfun.R")
#     and the rows correspond to the data groups.
nls.plot <- function(ds, c.df=NULL) {
  p <- ggplot(ds, aes(x=times, y=ft, col=grp)) +
    geom_point(size=2)
  if(is.null(c.df)) {
    p <- p + ggtitle("Synthetic data")
  } else {
    grpnames <- colnames(ds)
    for (ri in 1:nrow(c.df)) {
      p <- p + geom_function(fun=frap.fun, n=100, col="black",
                             args=as.list(c.df[ri, ]))
    }
    p <- p + ggtitle("FRAP curves fitted to the groups")
  }
  return(p)
}
