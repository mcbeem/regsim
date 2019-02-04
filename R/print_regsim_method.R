#' Print method for class regsim
#'
#' @param x an object of class regsim to be printed
#' @param ... additional arguments to be passed to \code{print}
#'
#' @export

print.regsim <- function (x, ...) {
  hid <- attr(x, "hidden")
  print(x[!names(x) %in% hid], ...)
}
