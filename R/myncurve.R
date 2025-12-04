#' myncurve
#'
#' @param mu Average Value
#' @param sigma Standard Deviation
#' @param a Interval Upper Limit
#'
#' @returns A list including entered values and probability, as well as a graph
#' @export
#'
#' @examples
#' myncurve(mu = 1, sigma = 2, a = 3)
myncurve = function(mu, sigma, a){
  curve(dnorm(x,mean=mu,sd=sigma), xlim = c(mu-3*sigma, mu + 3*sigma),col="Red", lwd=2,
        ylab="Normal Density", main=paste("Normal(", mu, ",", sigma, ")"))
  list(mu = mu, sigma = sigma)
  xcurve =  seq(mu-3*sigma, a, length = 1000)
  ycurve = dnorm(xcurve, mean=mu, sd=sigma)
  polygon(c(mu - 3*sigma, xcurve, a), c(0, ycurve, 0), col="Red")
  prob = pnorm(a, mean=mu, sd=sigma)
  prob = round(prob, 4)
return(list(mu = mu, sigma = sigma, a = a, prob = prob))
}
