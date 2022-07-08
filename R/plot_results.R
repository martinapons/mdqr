#' Plot results
#' @export
#'
#' @param object An object estimated by \code{mdqr}.
#' @param x A string containing the name of the regressor of interest
#' @return A \code{ggplot}
#' @import ggplot2
#' @section Details Additional ggplot option can be added after the command. For example: plot_results(md_object, "treatment") +  labs(x = "Quantiles", y = "Point Estimates" , title = "My Title")

plot_results <- function(object, x){
  tab <- result_table(object, x)
  tab <-  tab %>% dplyr::as_tibble() %>% dplyr::mutate(lb = Estimate - 1.96*`Std. Error`)
  tab <- tab %>% dplyr::mutate(ub = Estimate + 1.96*`Std. Error`)

  plot <-  ggplot2::ggplot(data = tab) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::aes(Quantiles, Estimate) +
    ggplot2::geom_line(size = 1, color = "dodgerblue4") +
    ggplot2::geom_point(color = "dodgerblue4", size = 2) +
    ggplot2::geom_ribbon(aes(ymin = lb, ymax = ub), alpha = .2, fill = "deepskyblue2")
  return(plot)
}
