#' Plot results
#' @export
#'
#' @param object An object of class \code{mdqr}.
#' @param x A string containing the name of the regressor of interest
#' @param level set confidence level; default is 95 percent.
#' @return A \code{ggplot}
#' @import ggplot2
#' @details # Additional options
#' Additional \code{\link[ggplot2]{ggplot}} options can be added after the command.
#' For example: \code{plot_results(md_object, "treatment") +  labs(x = "Quantiles", y = "Point Estimates" , title = "My Title")}


plot_mdqr <- function(object, x, level = 95){
  q <- 1- (1- 0.01*level)/2
  crit <- stats::qnorm(q)
  tab <- summary_mdqr(object, x)
  tab <-  tab %>% dplyr::as_tibble() %>% dplyr::mutate(lb = Estimate - crit*`Std. Error`)
  tab <- tab %>% dplyr::mutate(ub = Estimate + crit*`Std. Error`)

  plot <-  ggplot2::ggplot(data = tab) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::aes(Quantiles, Estimate) +
    ggplot2::geom_line(size = 1, color = "dodgerblue4") +
    ggplot2::geom_point(color = "dodgerblue4", size = 2) +
    ggplot2::geom_ribbon(aes(ymin = lb, ymax = ub), alpha = .2, fill = "deepskyblue2")
  return(plot)
}
