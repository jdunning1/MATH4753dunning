#' ntickets
#'
#' @param N The amount of seats on the plane
#' @param gamma The maximum allowed probability of overbooking the plane
#' @param p The probability that the person who bought a ticket shows up
#'
#' @returns A graph of the discrete binomial distribution, a graph of the normal
#' approximation and a list of nd, nc, N and gamma
#' @export
#'
#' @examples
#' ntickets(N=400,gamma = 0.02, p = 0.95)
ntickets <- function(N, gamma, p) {
  ### N = Number of seats on the plane
  ### gamma = Maximum probability that can be achieved of overbooking
  ### p = Probability that a passenger that bought a ticket shows up
  ### q = 1-p = Probability that a passenger that bought a ticket doesn't show up
  ### nd = Number of tickets to sell calculated using discrete binomial distribution
  ### nc = Number of tickets to sell calculated using normal approximation of the binomial
  ### X ~ bin(n, p)
  ### Start with binomial discrete
  ### Find the value of nd
  n <- N
  n1 <- N
  while(pbinom(N, n1, p) >= 1-gamma){
    n1 <- n1+1
  }
  nd <- n1-1

  ### Find the value of nc
  n2 <- N
  mu1 <- n2*p
  sigma1 <- sqrt(n2*p*(1-p))
  z <- (N+0.5-mu1)/sigma1
  while(pnorm(z) >= 1-gamma){
    mu1 <- n2*p
    sigma1 <- sqrt(n2*p*(1-p))
    z <- (N+0.5-mu1)/sigma1
    n2 <- n2+1
  }
  nc <- n2-1

  ### Find the Objective functions for both distributions
  min_val <- N
  max_val <- max(nd, nc)
  n_range <- seq(min_val, max_val + 30, by = 2)
  mu <- n_range*p
  sigma <- sqrt(n_range*p*(1-p))
  z2 <- (N + 0.5 - (n_range * p)) / sqrt(n_range * p * (1 - p))
  obj1 <- 1-gamma-pbinom(N, n_range, p)
  obj2 <- 1-gamma-pnorm(z2)
  par(mfrow = c(2, 1))

  ### Plot the first graph using the discrete binomial distribution
  title1 <- paste("Objective Vs n to find optimal tickets sold\n(",
    nd, ") gamma= ", gamma, " N= ", N, " discrete")
  plot1 <- plot(n_range, obj1, type = "p", col = "blue", pch = 16, main = title1,,
       xlab = "n", ylab = "Objective", xlim = c(N, max(nd, nc) + 30), ylim = c(0,1))

  ### Add lines in between the points
  lines(n_range, obj1, type = "l", col = "black", lwd = 1)

  ### Create the vertical and horizontal lines that pass through the y intercept
  abline(h = 0, col = "red", lwd = 2)
  abline(v = nd, col = "red", lwd = 2)

  ### Plot the second graph using the normal approximation
  title2 <- paste("Objective Vs n to find optimal tickets sold\n(",
                  nc, ") gamma= ", gamma, " N= ", N, " continuous")
  plot2 <- plot(n_range, obj2, type = "l", main = title2, xlab = "n", ylab = "Objective",
                xlim = c(N, max(nd, nc) + 30), ylim = c(0,1))
  ### Create the vertical and horizontal lines that pass through the y intercept
  abline(h = 0, col = "blue", lwd = 2)
  abline(v = nd, col = "blue", lwd = 2)

  ### Return the list after creating both plots
  rlist <- list(nd = nd, nc = nc, N = N, gamma = gamma)
  return(rlist)
}
