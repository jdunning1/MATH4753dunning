#' @title A Normal curve with area
#'
#' @param x The upper bound of the area
#' @param mu the mean
#' @param sigma the standard deviation
#'
#' @return A curve, shaded region and area
#' @export
#'
#' @examples
#' myncurve(x=5,mu =10, sigma =4)
myncurve = function(x,mu, sigma){
  curve(dnorm(x,mean=mu,sd=sigma),
        xlim   =   c(mu-3*sigma,   mu   + 3*sigma),
        ylab = "Normal density")

  xcurve = seq(mu-4*sigma,x, length = 1000)
  ycurve = dnorm(xcurve, mean = mu, sd = sigma)
  polygon(c(xcurve[1], xcurve, xcurve[1000]), c(0,ycurve,0),col = "Pink")
  area = pnorm(x, mu,sigma)
  area = round(area,4)
  list(area = area)
  }
