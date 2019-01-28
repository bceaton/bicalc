## ----setup, include=FALSE------------------------------------------------
library(bicalc)
knitr::opts_chunk$set(echo = TRUE)

## ---- random_diameters, echo=TRUE----------------------------------------
diameters = 2^rnorm(200, mean = 5, sd = 1)   #generate some log-normally distributed random numbers
head(diameters)   #print the first few data points to the console

## ---- wrangle_data, echo=TRUE--------------------------------------------
  test.cfd = MakeCFD(diameters)
  head(test.cfd)   #print the CFD to the console

## ------------------------------------------------------------------------
D50.estimate = WolmanCI(cfd = test.cfd, n = 200, probs = c(0.5))
D50.estimate   #print the results for D50 to the console


## ------------------------------------------------------------------------
D.estimates = WolmanCI(cfd = test.cfd, n = 200, probs = c(0.25, 0.5, 0.75))
D.estimates   #print the results for D25, D50 and D75 to the console


## ---- fig.width=5, fig.height=5------------------------------------------
  
#  create the standard plot using base graphics
  plot(test.cfd,  
     type ="b",
     pch = 21,
     log = "x",
     col = rgb(0,0,1,0.35),
     ylim = c(0, 1),
     xlab = "grain size (mm)",
     ylab = "Proportion Finer"
     )

#  create the confidence interval polygon
test.poly = PolyCI(cfd = test.cfd, n = 200, alpha = 0.05)

#  add it to the basic plot
polygon(test.poly, col = rgb(0,0,1,0.2), lty = 0)

## ------------------------------------------------------------------------


