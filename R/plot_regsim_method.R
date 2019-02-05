#' Plot method for class regsim
#'
#' @param x an object of class regsim to be printed
#' @param sizeMan numeric; controls scale of the variables and boxes, defaults to 10
#' @param edge.label.cex numeric; controls scale of the path coefficient. defaults to 1.5
#' @param fixedStyle numeric or vector of an lty and a color for the arrows. defaults to 1 for black, solid lines
#' @param ... additional arguments to be passed to \code{semPaths}
#'
#' @importFrom semPlot semPlotModel semPaths
#' @export

plot.regsim <- function(x, sizeMan=10, edge.label.cex=1.5,
                        fixedStyle=1, ...) {

  semPlot::semPaths(semPlot::semPlotModel(x$true.model), whatLabels="est",
           residuals=FALSE, sizeMan=sizeMan, edge.label.cex=edge.label.cex,
           fixedStyle=fixedStyle, ...)
}
