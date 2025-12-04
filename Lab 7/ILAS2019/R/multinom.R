#' @title Multinomial simulation
#'
#'@description This is a function from Lab 5
#' @param iter number of iterations
#' @param n number of trials
#' @param p probability vector
#'
#' @return A barplot and a table of relative frequencies
#' @export
#'
#' @examples
#' mymult()
mymult=function(iter=100,n=10, p=c(1,1,1,1)/4){
  # make a matrix to hold the samples
  #initially filled with NA's
  sam.mat=matrix(NA,nr=n,nc=iter, byrow=TRUE)
  #The number of categories is k
  k=length(p)
  # Make a matrix that will hold the frequencies in each sample
  tab.mat=matrix(NA,nr=k,nc=iter, byrow=TRUE)


  for(i in 1:iter){
    #Fill each column with a new sample
    sam.mat[,i]=sample(1:k,n,replace=TRUE, prob=p)
    #Collect all the frequencies of each of the k values
    tab.mat[,i]=table(factor(sam.mat[,i],levels=1:k))
  }
  # sum the frequecies
  freq=apply(tab.mat,1,sum)
  # put names to them
  names(freq)=1:k
  #create a barplot of refative freq
  barplot(freq/(n*iter),col=rainbow(k) )
  tab.mat
}
