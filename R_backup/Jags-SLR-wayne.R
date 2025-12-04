# Jags-ExampleScript.R
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis:
#   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
library(rjags)
install.packages(rjags)
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
#rm(list=ls())  # Careful! This clears all of R's memory!

# Load the functions used below:
source("DBDA2E-utilities.R") # Must be in R's current working directory.
require(rjags)               # Must have previously installed package rjags.

fileNameRoot="Jags-ExampleScript" # For output file names.

# Load the data:
#myData = read.csv("z15N50.csv") # Read data file; must be in curr. work. dir.
#y = myData$y        # The y values are in the column named y.
#y <- 4
#Ntotal = length(y)  # Compute the total number of flips.
spruce <- read.csv("SPRUCE.csv")
dataList = list(    # Put the information into a list.
  x = spruce$BHDiameter ,
  y = spruce$Height,
  n = 36
)

# Define the model:
modelString = "
model {
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)  # Response variable follows a normal distribution
    mu[i] <- beta0 + beta1 * x[i]  # Linear predictor
  }

  # Priors
  beta0 ~ dnorm(0, 0.001)  # Prior for intercept (weakly informative)
  beta1 ~ dnorm(0, 0.001)  # Prior for slope (weakly informative)
  sigma ~ dunif(0, 100)    # Prior for standard deviation (uniform on 0 to 100)
  tau <- 1 / (sigma * sigma)  # Precision (1/variance)

  # Derived quantity (optional)
  residual_variance <- 1 / tau
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Initialize the chains based on MLE of data.
# Option: Use single initial value for all chains:
#  thetaInit = sum(y)/length(y)
#  initsList = list( theta=thetaInit )
# Option: Use function that generates random values for each chain:


# Run the chains:
jagsModel = jags.model( file="TEMPmodel.txt" ,
                        n.chains=3 , n.adapt=500 ,data=dataList)
update( jagsModel , n.iter=500 )
codaSamples = coda.samples( jagsModel , variable.names=c("beta0", "beta1", "sigma") ,
                            n.iter=3334 )

save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )

# Examine the chains:
# Convergence diagnostics:
diagMCMC( codaObject=codaSamples , parName="beta1" )
saveGraph( file=paste0(fileNameRoot,"beta1Diag") , type="jpg" )
# Posterior descriptives:
openGraph(height=3,width=4)
par( mar = c(3.5,0.5,2.5,0.5) , mgp = c(2.25,0.7,0) )
plotPost( codaSamples[,"beta1"] , main="beta1" , xlab=bquote(beta[1]) )
saveGraph( file=paste0(fileNameRoot,"beta1Post") , type="bmp" )
# Re-plot with different annotations:
plotPost( codaSamples[,"beta1"] , main="Posterior beta1 by Wayne" , xlab=bquote(beta[1]) ,
          cenTend="mean" , compVal=0.5 , ROPE=c(5/12-0.02, 5/12+0.02) , credMass=0.95 )
saveGraph( file=paste0(fileNameRoot,"beta1") , type="bmp" )
