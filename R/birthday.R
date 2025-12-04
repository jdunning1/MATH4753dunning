#' Birthday problem probability
#'
#' Calculates the probability for the birthday problem.
#' @param x The number of people.
#' @return The probability of least two people sharing a birthday.
#' @export
#' @examples
#' birthday(23)
birthday <- function(x) {
  1 - exp(lchoose(365, x) + lfactorial(x) - x * log(365))
}
