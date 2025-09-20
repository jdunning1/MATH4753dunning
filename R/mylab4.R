#' Piecewise regression fit
#'
#' Fits a piecewise linear regression with one breakpoint.
#'
#' @param data A data frame containing x and y variables.
#' @param x Name of the predictor variable (string).
#' @param y Name of the response variable (string).
#' @param knot Knot value
#'
#' @return A fitted lm object.
#' @examples
#' data(spruce, package = "Intro2R")
#' spruce.df <- spruce
#' fit <- piecewise_reg(spruce.df, "BHDiameter", "Height", knot = 18)
#' summary(fit)
#' @export
piecewise_reg <- function(data, x, y, knot) {
  xvals <- data[[x]]
  yvals <- data[[y]]
  piece <- ifelse(xvals > knot, xvals - knot, 0)
  df <- data.frame(y = yvals, x = xvals, piece = piece)
  lm(y ~ x + piece, data = df)
}
library(s20x)
getwd()
data(spruce)
layout(matrix(1:3, nrow = 3, ncol = 1), heights = c(3, 3, 3))
trendscatter(x = spruce.df$Height, y = spruce.df$BHDiameter, f=0.5)
trendscatter(x = spruce.df$Height, y = spruce.df$BHDiameter, f=0.6)
trendscatter(x = spruce.df$Height, y = spruce.df$BHDiameter, f=0.7)
