#' CLT function
#'
#' @description This function creates a samlpling distribution of a statistic
#'
#' @param n sample size
#' @param iter number of iterations (bigger the better - choice depends on computing power and how long you wish to wait)
#' @param a first parameter of uniform
#' @param b second parameter of uniform
#'
#' @return Histogram and a vector of stats
#' @export
#'
#' @examples
#' myclt(n=10,iter=10000)
myclt=function(n,iter,a=0,b=1){
  y=runif(n*iter,a,b)
  data=matrix(y,nr=n,nc=iter,byrow=TRUE)
  sm=apply(data,2,sum)
  h=hist(sm,plot=FALSE)
  hist(sm,col=rainbow(length(h$mids)),freq=FALSE,main="Distribution of the sum of uniforms")
  curve(dnorm(x,mean=n*(a+b)/2,sd=sqrt(n*(b-a)^2/12)),add=TRUE,lwd=2,col="Blue")
  sm
}
