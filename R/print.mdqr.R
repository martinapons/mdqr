#' Print results of mdqr object
#' @param x An object of class mdqr.
#' @param ... Additional arguments to be passed to the print method.
#' @return The function prints the first element of the mdqr object and returns the object invisibly.
#' @export

print.mdqr <- function(x, ...){
print(x[[1]])
invisible(x)
}

