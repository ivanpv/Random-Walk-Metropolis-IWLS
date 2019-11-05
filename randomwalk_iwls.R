#RANDOM WALK M-H############################################################################################################

theta <- 6
count<-0
rwmet <- function(M=10000,phi=-2){
  #random walk function will let you iterate by entering a Sample size and a Phi value which will iterate through our proposal.
  vec<-vector("numeric", M)
  x<-1
  #for the proposal we are using the normal distribution but since the rnorm funcion is using the standard deviation we need to use square root on our exp(phi)
  proposed<-rnorm(1,x,sqrt(exp(phi)))
  acceptProb<-min(1, exp(theta*proposed-exp(theta*proposed) - (theta*x - exp(theta*x))))
  #notice the next conditional for initial value is to know if proposal passed the acceptance probability condition and add 1 to counter if so.
  if (acceptProb >= runif(1)){
    vec[1]<-proposed
    count<-count+1 
  }else{ 
    vec[1]<- x
  }
  for (i in 2:M) {
    #for the proposal we are using the normal distribution but since the rnorm funcion is using the standard deviation we need to use square root on our exp(phi)
    proposed<-rnorm(1,vec[i-1],sqrt(exp(phi)))
    acceptProb<-min(1, exp(theta*proposed-exp(theta*proposed) - (theta*vec[i-1] - exp(theta*vec[i-1]))))
    if (acceptProb >= runif(1)){
      count<-count+1
      vec[i]<-proposed
    }else{
      vec[i]<- vec[i-1]
    }
  }
  rate<-count/M
  list(vec,rate)
}
#To consider, there is a set.seed of 100 in order to have the same results everytime I run it.
set.seed(100)
print(rwmet(10000,-2)[1])


#Since we want to ensure a good approximation we are looking for a low variance at which we will
#obtained by looping through a 100 cases and saving their mean in a vector that I am creating and finally obtaining the variance of all of them.
result<-vector("numeric", 100)
set.seed(100)
for (N in 1:100) {
  result[N] <- mean((rwmet(10000,-4)[[1]]))
}
print(var(result))

set.seed(100)
M <-10000
C<-100
phiC <-seq(-5,10,length = C)
acceptRate <- rep(0,C)
#This loop will fill our acceptRate vector with the acceptance rate of their respective values from the range.
for (index in 1:C) {
  acceptRate[index]<-rwmet(M,phiC[index])[[2]]
}
plot(phiC,acceptRate)
axis(1,c(-5:10))
abline(h=.234, col="red")
legend(1.5, 0.5, legend = c("Optimal value appears", "to be 0 as shown"))


#Based on the plot we are expecting it to be 0 or very close to 0.s
count<-1
for (indexCount in acceptRate) {
  if (indexCount < .234){
      print(phiC[count])
      break()
  }else{
    count<-count+1
  }
}

#Finally we print out the line of our optimal value for phi which is 0.
abline(v=0, col="red")
#Looking at the plot the phi where the acceptance rate is closer to .234 is close to 0. By looping the index I find out that indeed it was 0 so the optimal value for phi will be 0
#for our theta = 6


#Now for this case we need to do the same thing as we did for Question 1 Exercise C. We will loop through our means in order to obtain a variance. But using our optimal phi value which is 0
result<-vector("numeric", 100)
set.seed(100)
for (N in 1:100) {
  result[N] <- mean(rwmet(10000,0)[[1]])
}
print(var(result))


#IWLS ######################################################################################################################

counter<-0
#Iterated weighted least squares function, will take binary vector and a model matrix to return the maximum likelihood estimates of Beta.
iwls<-function(X,y,itmax=25,eps=.00000001){
  size<-nrow(X)[1]
  W<-diag(0,size,size)
  z<-rep(0,size)
  #we will need to convert our information from coefficients of the file into a vector to be able to start the multiplication for mui that is why
  #we are using as.vector and lsfit in order to get least squares estimate of the B model of the input values to create our vector.
  b0<-as.vector(lsfit(X[,2:dim(X)[2]],y)$coefficients)
  #I will use a while so it iterates until the maximum number of iterations occurs at which it will iterate from calculating first our value of miu
  #then gMiu in order to obtain our W Matrix and finally our Z vector for the lenght in our row matrix. Which I will use to obtain the new beta.
  while (counter<itmax) {
    for (index in 1:size) {
      m<- exp(t(X[index,])%*%b0)/(1+exp(t(X[index,])%*%b0))
      gm<- 1/(m * (1-m))
      W[index,index]<- 1/(m* (1-m) * gm * gm)
      z[index]<- (y[index]-m) * (gm)
    }
    #we will keep counting until reaching max iteration
    counter<-counter+1
    betanew<-b0+solve( (t(X) %*% W %*% X)) %*% t(X) %*% W %*% z
    b0<-betanew
  }
  return(betanew)
}



diab <- as.matrix(read.table("diab.txt", header = TRUE))
X<-cbind(rep(1,nrow(diab)[1]),diab[,1:8])
y<-diab[,9]


#Running our function using our model matrix and binary vector respectively.
iwls(X,y)


