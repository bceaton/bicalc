#' Calculate the percentiles and confidence intervals
#'
#' \code{QuantCI} is a function that uses cumulative frequency distribution
#' data for Wolman or grid-by-number samples of bed surface texture to
#' estimate the value of user-specified percentiles. The function also uses
#' binomial theory to estimate the confidence intervals corresponding to a
#' user-specified confidence level, based on the number of grain size
#' measurements that were used to construct the cumulative frequency
#' distribution.
#'
#' The function returns a data frame listing the estimate of each quantile,
#' the upper limit, and the lower limit. This method is applicable to all
#' quantiles between the 5th and 95th percentile of the distirbution, but
#' not beyond that.
#'
#' @param cfd data frame providing a list grain sizes in the first variable
#' and the corresponding cumulative proportion finer in the second. The
#' grain sizes should be recorded in mm, and the proportion finer in [0,1].
#' @param n total number of observations upon which the cumulative
#' frequency distribution in \code{cfd} is based
#' @param probs numeric vector of probabilities with values in [0,1].
#' @param Np if the cumulative frequency data was a sub-sample from a finite
#' population of size Np, then the confidence range will be smaller and a
#' finite population correction is applied.
#' @param alpha  the desired confidence level for which to calculate a
#' confidence interval in [0,1].
#' @export
QuantCI = function(cfd, n, probs = seq(0.05,0.95,0.01), Np = NA, alpha = 0.05){
  tcrit = qt(1 - alpha/2, n)  #calculate t statistic
  #  calculate binomial standard deviation
  if (is.na(Np)){
    sigma = sqrt(n*probs*(1-probs))
  }else if(!is.na(Np)){
    fpc = sqrt ((Np - n) / (Np - 1))
    sigma = sqrt(n*probs*(1-probs))*fpc
  }
  # estimate percentiles
  phi = log2(cfd[[1]])
  X = cfd[[2]]
  estimate = 2^approx(x = X, y = phi, xout = probs)[[2]]
  p.upper = (n*probs + tcrit*sigma)/n
  upper = 2^approx(x = X, y = phi, xout = p.upper)[[2]]
  p.lower = (n*probs - tcrit*sigma)/n
  lower = 2^approx(x = X, y = phi, xout = p.lower)[[2]]

  results = data.frame(probs, estimate, lower, upper)

  return(results)
}


#' Generate a data frame containing cumulative frequency distribution data
#'
#' \code{MakeCFD} is a function that creates a data frame containing the
#' cumulative grain size distribution information for a sample of n observations
#' of b axis diameter (in mm).
#'
#' The function returns a data frame listing upper bounds of each size class
#' containing data, as well as the cumulative proportion of the observations
#' that fall below the given grain size. The output of this function is the
#' required input for other functions in this package. In order to retain
#' information about the sample size, there is the option to output the total
#' counts for each size class, as well.
#'
#' @param obs vector containing b-axis measurements of gravel bed surface
#' @param increment an optional parameter to control the size of the grain size
#' classes used for the analysis (1.0 = 1.0 phi intervals, 0.5 = 1/2 phi
#' intervals, 0.25 = 1/4 phi intervals, etc.). Default value is 0.5.
#' @param count optional flag to record the number of observations falling in
#' each size class as an additional variable in the output data frame
#' @param plot optional flag to produce a graph of the resulting grain size
#' distribution
#' @export
MakeCFD = function(obs, increment = 0.5, count = FALSE, plot = FALSE){
  n = length(obs)
  max.phi = ceiling(log2(max(obs)))
  min.phi = floor(log2(min(obs)))
  size = 2^seq(from = min.phi, to = max.phi, by = increment)
  results = hist(obs,
                 breaks = size,
                 plot = FALSE)
  cfd = data.frame(size)
  cfd$probs = c(0, cumsum(results$counts))/sum(results$counts)

  if(count == TRUE){
    cfd$counts = c(0,results$counts)
    cfd = (cfd[which(!cfd$counts == 0),])
  }

  if(plot == TRUE){
    plot(cfd[[1]], cfd[[2]],
         type = "b",
         col = "blue",
         log = "x",
         xlab = "grain size (mm)",
         ylab = "cum. prop. finer")
  }

  return(cfd)
}

#' Generate a polygon representing the confidence bounds for a grain size distribution
#'
#' \code{PolyCI} is a function that generates a data frame containing the
#' coordinates defining a polygon representing the confidence bounds around
#' a bed surface grain size distribution. The function contains an option to
#' generate a plot of the distribution showing the estimates of the percentiles
#' between the D5 and the D95, as well as the confidence limits about those
#' estimates.
#'
#' @param cfd data frame providing a list grain sizes in the first variable
#' and the corresponding cumulative proportion finer in the second. The
#' grain sizes should be recorded in mm, and the proportion finer in [0,1].
#' @param n total number of observations upon which the cumulative
#' frequency distribution in \code{cfd} is based
#' @param probs numeric vector of probabilities defining the polygon vertices,
#'  with values in [0,1].
#' @param Np if the cumulative frequency data was a sub-sample from a finite
#' population of size Np, then the confidence range will be smaller and a
#' finite population correction is applied.
#' @param alpha  the desired confidence level for which to calculate a
#' confidence interval in [0,1].
#' @param plot optional flag to produce a graph of the resulting grain size
#' distribution
#' @export
PolyCI = function(cfd, n, probs = seq(0.05, 0.95, 0.01), Np= NA, alpha = 0.05, plot = FALSE){
  tcrit = qt(1 - alpha/2, n)  #calculate t statistic
  #  calculate binomial standard deviation
  if (is.na(Np)){
    sigma = sqrt(n*probs*(1-probs))
  }else if(!is.na(Np)){
    fpc = sqrt ((Np - n) / (Np - 1))
    sigma = sqrt(n*probs*(1-probs))*fpc
  }
  # estimate percentiles
  phi = log2(cfd[[1]])
  X = cfd[[2]]
  estimate = 2^approx(x = X, y = phi, xout = probs)[[2]]
  p.upper = (n*probs + tcrit*sigma)/n
  upper = 2^approx(x = X, y = phi, xout = p.upper)[[2]]
  p.lower = (n*probs - tcrit*sigma)/n
  lower = 2^approx(x = X, y = phi, xout = p.lower)[[2]]

  x.poly = c(upper, rev(lower))
  y.poly = c(probs, rev(probs))
  poly.out = data.frame(x.poly, y.poly)

  if(plot == TRUE){
    plot(cfd[[1]], cfd[[2]],
         type = "b",
         pch = 20,
         col = "blue",
         log = "x",
         xlim = c(min(cfd[[1]]),max(cfd[[1]])),
         ylim = c(0.05,0.95),
         xlab = "grain size (mm)",
         ylab = "cum. prop. finer")
    polygon(poly.out,
            col=rgb(0, 0, 1,0.3),
            lty = 0)
    abline(h = c(0.05, 0.95))
  }
  return(poly.out)
}
