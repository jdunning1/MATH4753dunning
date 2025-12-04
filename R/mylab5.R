## Make a bin to store each iteration of the matrix and then create a boxplot using said bin
#' mylab5
#'
#' @param iter # How many times to iterate
#' @param n # Number of values
#' @param p # Probability
#'
#' @return # A table and plot
#' @export
#'
#' @examples
#' mylab5(100, 15, 0.5)
mybin <- function(iter = 100, n = 10, p = 0.7) {
  sam.mat <- matrix(NA, nr = n, nc = iter, byrow = TRUE)
  succ <- c()
  for (i in 1:iter) {
    sam.mat[, i] <- sample(c(1, 0), n, replace = TRUE, prob = c(p, 1 - p))
    succ[i] <- sum(sam.mat[, i])
  }
  succ.tab <- table(factor(succ, levels = 0:n))
  barplot(succ.tab / iter, col = rainbow(n + 1), main = "Binomial simulation",
          xlab = "Number of successes", ylab = "Proportion")
  succ.tab / iter
}
