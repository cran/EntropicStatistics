#' Heat Map for Distribution Visualization
#'
#' Return a heat map that displays distributional characteristics for selected groups.
#'
#' @param data_frequency_list A \code{list} of nonnegative frequency/count vectors
#'   (e.g., from \code{table()}). Element names (if present) are used as row labels.
#' @param orders Numeric vector of generalized Shannon entropy orders to evaluate.
#' @param selection Integer indices selecting which elements of \code{data_frequency_list} to include.
#' @param plot_order Integer indices giving the order (bottom to top) of the selected groups.
#' @param RowNames Character vector of row labels for the selected groups
#'   (default: \code{names(data_frequency_list)[plot_order]}).
#' @param title Character string, plot title.
#' @param x_ticks Numeric vector of x-axis tick locations; values must be present in \code{orders}.
#' @param plot_margin Plot margin, a \code{ggplot2::margin()} object.
#' @param text_face Integer font face in the plot: \code{1} = "plain", \code{2} = "italic",
#'   \code{3} = "bold", \code{4} = "bold.italic".
#' @param fill_colors Character vector of three colors for low, mid, and high values.
#' @param title_text_size Numeric size of the title text.
#' @param label_text_size Numeric size of axis text.
#'
#' @details
#' Provides a quick, nonparametric view of distributional differences across multiple
#' groups simultaneously, using generalized Shannon entropy over a range of orders.
#' Input vectors should be nonnegative counts. Zero-probability categories are handled internally.
#'
#' @return A \code{ggplot} object representing the heat map.
#'
#' @references
#' Zhang, J. and Shi, J. (2024). Nonparametric clustering of discrete probability
#' distributions with generalized Shannon's entropy and heatmap. \emph{Statistics & Probability Letters}.
#' \doi{10.1016/j.spl.2024.110070}
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_tile}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' binom_n <- 10
#' sample_size <- 400
#' sample_1 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.1))
#' sample_2 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.2))
#' sample_3 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.3))
#' sample_4 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.4))
#' sample_5 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.5))
#' sample_6 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.6))
#' sample_7 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.7))
#' sample_8 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.8))
#' sample_9 <- table(stats::rbinom(n = sample_size, size = binom_n, prob = 0.9))
#' poisson_1 <- table(stats::rpois(n = sample_size, lambda = 1))
#' poisson_2 <- table(stats::rpois(n = sample_size, lambda = 2))
#' poisson_3 <- table(stats::rpois(n = sample_size, lambda = 3))
#' poisson_4 <- table(stats::rpois(n = sample_size, lambda = 4))
#' poisson_5 <- table(stats::rpois(n = sample_size, lambda = 5))
#' poisson_6 <- table(stats::rpois(n = sample_size, lambda = 6))
#' poisson_7 <- table(stats::rpois(n = sample_size, lambda = 7))
#' poisson_8 <- table(stats::rpois(n = sample_size, lambda = 8))
#' poisson_9 <- table(stats::rpois(n = sample_size, lambda = 9))
#' data_samples <- list(
#'   binom_0.1 = sample_1, binom_0.2 = sample_2, binom_0.3 = sample_3,
#'   binom_0.4 = sample_4, binom_0.5 = sample_5, binom_0.6 = sample_6,
#'   binom_0.7 = sample_7, binom_0.8 = sample_8, binom_0.9 = sample_9,
#'   Poisson_1 = poisson_1, Poisson_2 = poisson_2, Poisson_3 = poisson_3,
#'   Poisson_4 = poisson_4, Poisson_5 = poisson_5, Poisson_6 = poisson_6,
#'   Poisson_7 = poisson_7, Poisson_8 = poisson_8, Poisson_9 = poisson_9
#' )
#' HeatMap(data_samples)
#' HeatMap(data_samples, selection = sample(seq_along(data_samples), 6))
#' HeatMap(data_samples, selection = 1:9)
#' HeatMap(data_samples, selection = 10:13)
#' HeatMap(data_samples, selection = 14:18)
#' }
#'
#' @concept heatmap
#' @concept entropy
#' @keywords hplot
#'
#' @importFrom ggplot2 ggplot aes element_text margin geom_tile scale_fill_gradientn
#' @importFrom ggplot2 scale_x_discrete scale_y_discrete ggtitle theme
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
HeatMap <- function(
    data_frequency_list,
    orders = seq(0.50, 3, by = 0.01),
    selection = 1:length(data_frequency_list),
    plot_order = selection,
    RowNames = names(data_frequency_list)[plot_order],
    title = "HeatMap",
    x_ticks = round(stats::quantile(orders, c(0, 0.25, 0.5, 0.75, 1)), 2),
    plot_margin = ggplot2::margin(0.5, 0.2, 0.2, 1, "cm"),
    text_face = 1,
    fill_colors = c("blue4", "white", "red3"),
    title_text_size = 25,
    label_text_size = 25
) {
  generalized.entropy <- function(sample_freq, order) {
    n <- sum(sample_freq)
    if (n <= 0) return(rep(NA_real_, length(order)))
    phats <- sample_freq / n
    phats <- phats[phats > 0]
    vapply(order, function(o) {
      pmhats <- phats^o / sum(phats^o)
      sum(-pmhats * log(pmhats))
    }, numeric(1))
  }

  df <- t(data.frame(lapply(
    data_frequency_list[plot_order],
    function(x) generalized.entropy(x, orders)
  )))
  rownames(df) <- RowNames
  colnames(df) <- orders

  p <- as.data.frame(df) |>
    tibble::rownames_to_column(var = "rowname") |>
    tidyr::pivot_longer(cols = - "rowname", names_to = "name", values_to = "value") |>
    dplyr::mutate(rowname = factor(.data$rowname, unique(rownames(df)))) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = factor(.data$name, unique(.data$name)),
        y = .data$rowname,
        fill = .data$value
      )
    ) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_discrete(
      name = "Orders",
      breaks = as.character(x_ticks),
      labels = function(x) x
    ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = fill_colors) +
    ggplot2::scale_y_discrete(name = NULL, position = "right") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = title_text_size),
      text       = ggplot2::element_text(size = label_text_size, face = text_face),
      legend.position = "none",
      axis.ticks = ggplot2::element_blank(),
      plot.margin = plot_margin
    )

  return(p)
}
