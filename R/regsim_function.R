#' Fit and appraise regression models using simulated data
#'
#' This function provides a quick and convenient means for simulating data from a
#' path diagram or SEM, fitting a regression model to the simulated data, and appraising the
#' model's performance with respect to a particular target parameter and its standard
#' error. The function calculates the bias, coverage rate, RMSE, and analytical and empirical
#' standard errors.
#'
#' @param reps the number of repetitions to perform
#' @param n the sample size of the generated data
#' @param true.model character; the population model used to generate the data specified in \code{lavaan} syntax. See example.
#' @param fit.model character; the model to the fitted to the data in \code{lm()} model formula syntax
#' @param targetparm character; the focal predictor. Must be included in the model formula passed to \code{true.model}
#' @param targetval the true value of the path coefficient relating the focal predictor to the response variable
#' @param interval the confidence interval width, defaults to 0.95.
#' @param standardized logical; should the data be generated on a standardized scale? defaults to FALSE
#' @param ... additional arguments passed to lavaan's \code{simulateData} function. For example, one could invoke the
#' \code{kurtosis} or \code{skewness} arguments to generate non-normal data
#'
#' @return A list containing the following:
#'  \item{true.model}{the true model from which data are simulated}
#'  \item{fit.model}{the regression model that is fit to the data}
#'  \item{b}{a vector of target parameter estimates across the repetitions}
#'  \item{data}{a data frame containing the simulated data from the first repetition}
#'  \item{true.DAG}{the true.model expressed as an object of class dagitty}
#'  \item{targetval}{the true value of the target parameter, used to calculate bias, RMSE, and coverage}
#'  \item{expected.b}{the mean value of the target parameter across the repetitions}
#'  \item{bias}{the mean difference between the estimated and target values of the parameter across repetitions}
#'  \item{analytic.SE}{the mean of the analytic standard errors across repetitions}
#'  \item{empirical.SE}{the standard deviation of the parameter estimates}
#'  \item{analytic.CI}{the analytic confidence interval boundaries, calculated as the mean of the lower and upper boundaries across the repetitions}
#'  \item{empirical.CI}{the empirical confidence interval boundaries across the repetitions}
#'  \item{coverage}{the proportion of repetitions in which the confidence interval captured
#'  the specified target value of the parameter. This should appoximately equal
#'  the specified confidence level}
#'  \item{RMSE}{the root mean squared error}
#'  \item{adjustment.sets}{the sets of covariates that must be included in the fitted model
#'    to yield an unbiased and consistent estimate of the effect of the target parameter on the
#'    response variable in fit.model}
#'
#' @examples
#'
#' #  this is a DAG where Z is a confounder of the X -> Y relationship
#' true.model <- 'X ~ .5*Z
#'                Y ~ .5*X + .5*Z'
#'
#' # misspecified model omitting Z
#' fit.model <- 'Y ~ X'
#'
#' result <- regsim(reps=1000, n=100, true.model, fit.model, targetparm="X",
#'                  targetval=.5)
#' result
#'
#' # visualize the model
#' plot(result, type="model")
#' # visualize the simulated data
#' plot(result, type="data")
#' # visualize the results
#' plot(result, type="perf")
#' plot(result, type="compareSE")
#'
#' @importFrom stats as.formula confint lm na.omit sd quantile density dnorm
#' @importFrom lavaan simulateData
#' @importFrom graphics abline hist legend plot points
#' @export

regsim <- function(reps, n, true.model, fit.model, targetparm, targetval,
                   interval=.95, standardized=FALSE, ...) {

  # begin checks of the supplied arguments
    # check whether reps, n, and interval are numeric
    if (!is.numeric(n)) stop("n must be a number")
    if (!is.numeric(reps)) stop("reps must be a number")
    if (!is.numeric(interval)) stop("reps must be a number")

    # check whether standardized is logical
    if (!is.logical(standardized)) stop("standardized must be TRUE or FALSE")

    # check whether reps > 2
    if (reps < 2) stop("reps must be at least 2, and should be at least 1000 for trustworthy results")
    if (reps < 1000) warning("reps should be at least 1000 for trustworthy results")

    # check whether interval is a proportion
    if (0 > interval | interval > 1) stop("interval must be a proportion")

    # check if n is too small given fit.model
    if (n < length(all.vars(as.formula(fit.model)))+1) stop("insufficient df to estimate model; increase n")

    #check whether targetparm is a character
    if (is.character(targetparm) != TRUE) stop("The targetparm argument must be a quoted variable name")

    #check whether the targetparm is included in the fit.model and true.model
    if (!targetparm %in% all.vars(as.formula(fit.model))) stop("The target parameter must be included in fit.model")
    if (!grepl(targetparm, true.model, fixed=TRUE)) stop("The target parameter must be included in true.model")

    #check whether fit.model is a subset of true.model
    if (min(all.vars(as.formula(fit.model)) %in% unique(na.omit(unlist(strsplit(unlist(true.model), "[^a-zA-Z_.-]+")))))==0) stop("All variables in fit.model must be included in true.model")
  # end argument checks

  # find position of target parameter in model formula.
  parmnumber <- which(all.vars(as.formula(fit.model))==targetparm)

  # generate data using lavaan
  #  note that it is one huge dataset. this is much faster
  myData <- lavaan::simulateData(true.model, sample.nobs=n*reps,
                                 standardized=standardized, ...)

  # define empty objects to hold results
  coverage <- rep(NA, times=reps)
  bias <- rep(NA, times=reps)
  b <- rep(NA, times=reps)
  se <- rep(NA, times=reps)
  reg <- list(length=reps)
  CI <- matrix(c(NA, NA), nrow=reps, ncol=2)

  # loop over reps
  for (i in 1:reps) {
    # indexes subsets of myData for each rep
    reg[[i]] <- lm(fit.model, data=myData[(1+(i-1)*n):(i*n),])
    # b is the target parameter estimate
    b[i] <- reg[[i]]$coef[parmnumber]
    # se is the analytical standard error of the target estimate
    se[i] <- summary(reg[[i]])$coef[parmnumber, 2]
    # CI will contain the boundaries of the CI
    CI[i,] <- confint(reg[[i]], parm=parmnumber, level=interval)
    # coverage is a logical vector indicating whether true value was captured
    coverage[i] <- CI[i, 1] <= targetval & CI[i, 2] >= targetval
    # bias is the true value minus the estimate
    bias[i] <- reg[[i]]$coef[parmnumber] - targetval
  }

  true.DAG <- dagitty::lavaanToGraph(lavaan::lavaanify(true.model))
  dagitty::exposures(true.DAG) <- targetparm
  dagitty::outcomes(true.DAG) <- all.vars(as.formula(fit.model))[1]

  analytic.CI <- c(mean(CI[,1]), mean(CI[,2]))
  names(analytic.CI) <- c(paste0((1-interval)/2*100, "%"), paste0((1-((1-interval)/2))*100, "%"))

  output <- list(true.model=true.model,
                 fit.model=fit.model,
                 b=b,
                 data=myData[1:n,],
                 true.DAG = true.DAG,
                 targetval=targetval,
                 expected.b=mean(b),
                 bias=mean(bias),
                 analytic.SE=mean(se),
                 empirical.SE=sd(b),
                 analytic.CI=analytic.CI,
                 empirical.CI=quantile(b, c((1-interval)/2, 1-((1-interval)/2))),
                 coverage=mean(coverage),
                 RMSE = sqrt(mean(bias^2)),
                 adjustment.sets = dagitty::adjustmentSets(true.DAG))

  class(output) = c("regsim", "list")
  attr(output, "hidden") <- c("data", "b", "true.model", "fit.model", "true.DAG")


  # calculate the summary results
  return(output)
}
