#' Plot method for class regsim
#'
#' @param x an object of class regsim to be printed
#' @param type character; one of "model", "perf", or "data". the type of plot to produce.
#'     defaults to "model"
#' @param sizeMan numeric; controls scale of the variables and boxes, defaults to 10
#' @param edge.label.cex numeric; controls scale of the path coefficient. defaults to 1.5
#' @param fixedStyle numeric or vector of an lty and a color for the arrows. defaults to 1 for black, solid lines
#' @param nCharNodes numeric; controls how many characters of each variable's name to display. defaults to
#'    zero so that the full name is shown
#' @param ... additional arguments to be passed to \code{semPlot::semPaths} (for \code{type="model"}),
#'    \code{hist} (for \code{type="perf"}), or to \code{psych::pairs.panels} (for \code{type="data"})
#'
#' @importFrom semPlot semPlotModel semPaths
#' @export

plot.regsim <- function(x, type="model", sizeMan=10, edge.label.cex=1.5,
                        fixedStyle=1, nCharNodes=0, ...) {

  if (max(!type %in% c("model", "perf", "data", "compareSE"))) stop("Argument type must only contain \"model\", \"perf\", \"data\", or \"compareSE\"")

  if (max(type %in% "model")) {
    semPlot::semPaths(semPlot::semPlotModel(x$true.model), whatLabels="est",
            residuals=FALSE, exoCov=FALSE, sizeMan=sizeMan, edge.label.cex=edge.label.cex,
             fixedStyle=fixedStyle, ...)}

  if (max(type %in% "perf")) {
    lower.limit <- min(x$targetval, x$b)
    upper.limit <- max(x$targetval, x$b)
    hist(x$b, freq=FALSE, xlab="Estimates", main="",
         xlim=c(lower.limit, upper.limit), ...)
    points(density(x$b), type='l')
    abline(v=x$targetval, col="red", lty="dashed", lwd=2)
    abline(v=x$empirical.CI, lty="dotted")
  }

  if (max(type %in% "data")) {
    psych::pairs.panels(x$data, ...)
  }

  if (max(type %in% "compareSE")) {
    lower.limit <- min(x$targetval, x$b)
    upper.limit <- max(x$targetval, x$b)

    xvals <- seq(lower.limit, upper.limit, length.out=1000)
    yvals <- dnorm(xvals, mean=x$expected.b, sd=x$analytic.SE)
    dens <- density(x$b)
    plot(dens, ylim=c(0, max(dens$y, yvals)+.25*abs(max(dens$y, yvals))),
         xlab="Estimates", main="")
    points(xvals, yvals, type="l", lty="dotted", col="blue")
    abline(v=x$targetval, col="red", lty="dashed")
    legend("topright", legend=c("Empirical", "Analytic"),
           col=c("black", "blue"), lty=c("solid", "dashed"), cex=0.85)
  }
}
