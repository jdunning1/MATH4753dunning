# Shiny app: Demonstrating MLE for multiple univariate distributions
# Distributions included: Normal, Exponential, Poisson, Gamma, Beta
# Usage: save this file as app.R and run `shiny::runApp('.')` or open in RStudio and click Run App

library(shiny)
library(ggplot2)
library(MASS)     # fitdistr
library(stats4)   # optional

# Helper: negative log-likelihood for Beta parameterized by log(a), log(b)
negloglik_beta_trans <- function(logpar, x){
  a <- exp(logpar[1]); b <- exp(logpar[2])
  # ensure x in (0,1)
  if(any(x<=0) || any(x>=1)) return(1e12)
  -sum(dbeta(x, shape1=a, shape2=b, log=TRUE))
}

ui <- fluidPage(
  titlePanel("MLE playground: 5+ univariate distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Distribution",
                  choices = c("Normal", "Exponential", "Poisson", "Gamma", "Beta"),
                  selected = "Normal"),
      numericInput("n", "Sample size", 200, min = 10, max = 10000, step = 10),
      numericInput("seed", "Random seed (0 = random)", 0),
      uiOutput("param_inputs"),
      actionButton("gen", "(Re)generate data"),
      hr(),
      checkboxInput("show_ll_surface", "Show likelihood surface (if 2 params)", TRUE),
      helpText("This app simulates data from a chosen distribution, performs maximum likelihood estimation, and visualizes results.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("mainplot", height = "500px")),
        tabPanel("Likelihood", plotOutput("llplot", height = "500px")),
        tabPanel("Summary", verbatimTextOutput("summary")),
        tabPanel("Data (first rows)", tableOutput("head"))
      )
    )
  )
)

server <- function(input, output, session){
  # dynamic parameter inputs
  output$param_inputs <- renderUI({
    switch(input$dist,
           "Normal" = tagList(
             numericInput("norm_mean", "True mu", 0),
             numericInput("norm_sd", "True sigma", 1, min = 1e-6)
           ),
           "Exponential" = numericInput("exp_rate", "True rate (lambda)", 1, min = 1e-6),
           "Poisson" = numericInput("pois_lambda", "True lambda", 3, min = 0),
           "Gamma" = tagList(
             numericInput("gamma_shape", "True shape (k)", 2, min = 1e-6),
             numericInput("gamma_rate", "True rate (theta^-1)", 1, min = 1e-6)
           ),
           "Beta" = tagList(
             numericInput("beta_a", "True a (alpha)", 2, min = 1e-6),
             numericInput("beta_b", "True b (beta)", 5, min = 1e-6)
           )
    )
  })

  # reactive data
  dat <- eventReactive(input$gen, {
    if(input$seed != 0) set.seed(input$seed)
    n <- input$n
    switch(input$dist,
           "Normal" = rnorm(n, mean = input$norm_mean, sd = input$norm_sd),
           "Exponential" = rexp(n, rate = input$exp_rate),
           "Poisson" = rpois(n, lambda = input$pois_lambda),
           "Gamma" = rgamma(n, shape = input$gamma_shape, rate = input$gamma_rate),
           "Beta" = {
             # avoid exact 0/1 to keep log-liks finite
             x <- rbeta(n, shape1 = input$beta_a, shape2 = input$beta_b)
             pmin(pmax(x, 1e-6), 1-1e-6)
           }
    )
  }, ignoreNULL = FALSE)

  # MLE computation
  mle_res <- reactive({
    x <- dat()
    req(x)
    dist <- input$dist
    out <- list()

    if(dist == "Normal"){
      fit <- tryCatch(fitdistr(x, "normal"), error = function(e) NULL)
      out$method <- "fitdistr (normal)"
      out$coef <- if(!is.null(fit)) fit$estimate else c(mean=mean(x), sd=sd(x))
      out$se <- if(!is.null(fit)) fit$sd else c(NA, NA)
      out$loglik <- if(!is.null(fit)) fit$loglik else sum(dnorm(x, mean=out$coef[1], sd=out$coef[2], log=TRUE))
      out$npar <- 2
    } else if(dist == "Exponential"){
      fit <- tryCatch(fitdistr(x, "exponential"), error = function(e) NULL)
      out$method <- "fitdistr (exponential)"
      # fitdistr returns 'rate'
      out$coef <- if(!is.null(fit)) fit$estimate else c(rate=1/mean(x))
      out$se <- if(!is.null(fit)) fit$sd else c(NA)
      out$loglik <- if(!is.null(fit)) fit$loglik else sum(dexp(x, rate=out$coef[1], log=TRUE))
      out$npar <- 1
    } else if(dist == "Poisson"){
      fit <- tryCatch(fitdistr(x, "Poisson"), error = function(e) NULL)
      out$method <- "fitdistr (Poisson)"
      out$coef <- if(!is.null(fit)) fit$estimate else c(lambda=mean(x))
      out$se <- if(!is.null(fit)) fit$sd else c(NA)
      out$loglik <- if(!is.null(fit)) fit$loglik else sum(dpois(x, lambda=out$coef[1], log=TRUE))
      out$npar <- 1
    } else if(dist == "Gamma"){
      # fitdistr for gamma returns shape and rate
      fit <- tryCatch(fitdistr(x, "gamma"), error = function(e) NULL)
      out$method <- "fitdistr (gamma)"
      out$coef <- if(!is.null(fit)) fit$estimate else c(shape=NA, rate=NA)
      out$se <- if(!is.null(fit)) fit$sd else c(NA, NA)
      out$loglik <- if(!is.null(fit)) fit$loglik else sum(dgamma(x, shape=out$coef[1], rate=out$coef[2], log=TRUE))
      out$npar <- 2
    } else if(dist == "Beta"){
      # custom MLE: optimize over log(a), log(b)
      x2 <- pmin(pmax(x, 1e-6), 1-1e-6)
      start <- c(log(mean(x2)*( (mean(x2)*(1-mean(x2))) / var(x2) -1 )),
                 log((1-mean(x2))*( (mean(x2)*(1-mean(x2))) / var(x2) -1 )))
      start[is.na(start) | !is.finite(start)] <- log(c(1,1))
      opt <- optim(start, fn = negloglik_beta_trans, x = x2, hessian = TRUE, control = list(maxit=2000))
      logpar_hat <- opt$par
      a_hat <- exp(logpar_hat[1]); b_hat <- exp(logpar_hat[2])
      # delta method for se: var(a) ≈ a^2 * var(loga)
      cov_phi <- tryCatch(solve(opt$hessian), error=function(e) matrix(NA,2,2))
      se_a <- if(all(is.finite(cov_phi))) sqrt((a_hat^2) * cov_phi[1,1]) else NA
      se_b <- if(all(is.finite(cov_phi))) sqrt((b_hat^2) * cov_phi[2,2]) else NA
      out$method <- "optim on log-params (Beta)"
      out$coef <- c(a = a_hat, b = b_hat)
      out$se <- c(a = se_a, b = se_b)
      out$loglik <- -opt$value
      out$npar <- 2
      out$optim <- opt
    }

    out
  })

  output$mainplot <- renderPlot({
    x <- dat()
    req(x)
    d <- data.frame(x = x)
    dist <- input$dist
    est <- mle_res()

    p <- NULL
    if(dist %in% c("Normal", "Gamma", "Exponential", "Beta")){
      # continuous: histogram + fitted density
      p <- ggplot(d, aes(x)) + geom_histogram(aes(y=..density..), bins = 30, alpha=0.4) + ggtitle(paste0("Data and fitted density (", dist, ")"))
      if(dist == "Normal"){
        mu <- est$coef[1]; s <- est$coef[2]
        dens <- data.frame(x = seq(min(x), max(x), length=400))
        dens$y <- dnorm(dens$x, mean=mu, sd=s)
        p <- p + geom_line(data=dens, aes(x=x,y=y), size=1)
      } else if(dist == "Exponential"){
        rate <- est$coef[1]
        dens <- data.frame(x = seq(0, max(x), length=400))
        dens$y <- dexp(dens$x, rate=rate)
        p <- p + geom_line(data=dens, aes(x=x,y=y), size=1)
      } else if(dist == "Gamma"){
        shape <- est$coef[1]; rate <- est$coef[2]
        dens <- data.frame(x = seq(min(x), max(x), length=400))
        dens$y <- dgamma(dens$x, shape=shape, rate=rate)
        p <- p + geom_line(data=dens, aes(x=x,y=y), size=1)
      } else if(dist == "Beta"){
        a <- est$coef[1]; b <- est$coef[2]
        dens <- data.frame(x = seq(0,1,length=400))
        dens$y <- dbeta(dens$x, shape1=a, shape2=b)
        p <- p + geom_line(data=dens, aes(x=x,y=y), size=1)
        p <- p + xlim(0,1)
      }
    } else if(dist == "Poisson"){
      # discrete: barplot of frequencies + fitted Poisson PMF
      freq <- as.data.frame(table(x))
      freq$x <- as.numeric(as.character(freq$x))
      p <- ggplot(freq, aes(x=x, y=Freq/nrow(d))) + geom_col() + ggtitle("Empirical pmf and fitted Poisson") + ylab("Probability")
      lambda <- est$coef[1]
      xmax <- max(freq$x, qpois(0.999, lambda=lambda))
      pmf <- data.frame(x = 0:xmax)
      pmf$pmf <- dpois(pmf$x, lambda=lambda)
      p <- p + geom_point(data=pmf, aes(x=x,y=pmf), color = "black", size=2)
    }

    print(p)
  })

  output$llplot <- renderPlot({
    x <- dat(); req(x)
    est <- mle_res(); req(est)
    dist <- input$dist

    if(!input$show_ll_surface){
      plot.new(); text(0.5,0.5,"Likelihood surface hidden by user")
      return()
    }

    if(est$npar == 1){
      # 1-d log-likelihood curve
      if(dist == "Exponential"){
        rgrid <- seq(max(1e-6, est$coef[1]*0.2), est$coef[1]*3 + 1e-6, length=300)
        ll <- sapply(rgrid, function(rr) sum(dexp(x, rate=rr, log=TRUE)))
        plot(rgrid, ll, type='l', xlab='rate', ylab='log-likelihood', main='Log-likelihood (1D)')
        abline(v=est$coef[1], col='red')
      } else if(dist == "Poisson"){
        lgrid <- seq(max(1e-6, est$coef[1]*0.2), est$coef[1]*3 + 1e-6, length=300)
        ll <- sapply(lgrid, function(lam) sum(dpois(x, lambda=lam, log=TRUE)))
        plot(lgrid, ll, type='l', xlab='lambda', ylab='log-likelihood', main='Log-likelihood (1D)')
        abline(v=est$coef[1], col='red')
      }
    } else if(est$npar == 2){
      # 2-d likelihood surface (contour)
      if(dist == "Normal"){
        mu_hat <- est$coef[1]; s_hat <- est$coef[2]
        mu_grid <- seq(mu_hat - 2*sd(x), mu_hat + 2*sd(x), length=60)
        s_grid <- seq(max(1e-6, s_hat*0.2), s_hat*3, length=60)
        llmat <- outer(mu_grid, s_grid, Vectorize(function(mu,s) sum(dnorm(x, mean=mu, sd=s, log=TRUE))))
        contour(mu_grid, s_grid, llmat, xlab='mu', ylab='sigma', main='Log-likelihood contour (Normal)')
        points(mu_hat, s_hat, col='red', pch=19)
      } else if(dist == "Gamma"){
        sh_hat <- est$coef[1]; r_hat <- est$coef[2]
        sh_grid <- seq(max(1e-3, sh_hat*0.2), sh_hat*3, length=50)
        r_grid <- seq(max(1e-3, r_hat*0.2), r_hat*3, length=50)
        llmat <- outer(sh_grid, r_grid, Vectorize(function(sh,rr) sum(dgamma(x, shape=sh, rate=rr, log=TRUE))))
        contour(sh_grid, r_grid, llmat, xlab='shape', ylab='rate', main='Log-likelihood contour (Gamma)')
        points(sh_hat, r_hat, col='red', pch=19)
      } else if(dist == "Beta"){
        # use log-params for stability
        a_hat <- est$coef[1]; b_hat <- est$coef[2]
        a_grid <- exp(seq(log(max(1e-3, a_hat*0.2)), log(a_hat*3 + 1e-6), length=60))
        b_grid <- exp(seq(log(max(1e-3, b_hat*0.2)), log(b_hat*3 + 1e-6), length=60))
        llmat <- outer(a_grid, b_grid, Vectorize(function(a,b) sum(dbeta(x, shape1=a, shape2=b, log=TRUE))))
        contour(a_grid, b_grid, llmat, xlab='a', ylab='b', main='Log-likelihood contour (Beta)')
        points(a_hat, b_hat, col='red', pch=19)
      }
    }
  })

  output$summary <- renderPrint({
    est <- mle_res(); req(est)
    cat("Method:", est$method, "\n\n")
    cat("MLE estimates:\n")
    print(est$coef)
    cat("\nStandard errors (approx):\n")
    print(est$se)
    cat("\nLog-likelihood:", est$loglik, "\n")
    if(!is.null(est$optim)){
      cat("\nOptimizer output (Beta):\n")
      print(est$optim)
    }
  })

  output$head <- renderTable({
    head(dat())
  })
}

shinyApp(ui, server)
