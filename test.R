theta <- 5
xo <- 1
M <- 10000
xis<-rep(0,M)
xis[1]<- xo
phi <- 2

rwmet <- function(M,phi){
  theta <- 6
  xo <- 1
  M <- 10000
  xis<-rep(0,M)
  xis[1]<- xo
  for (i in 2:M) {
    y<-rnorm(1,xis[i-1],sqrt(exp(phi)))
    alph<-min(1,exp(theta*y-exp(theta*y)-(theta*xis[i-1]-exp(theta*xis[i-1]))))
    u<-runif(1)
    if(alph>=u){
      xis[i]<-y
    }
    else{
      xis[i]<-xis[i-1]
    }
  }
  store_vector<-xis
}
set.seed(1)
print(rwmet(M,phi))
print(var(rwmet(M,phi)))


#######
C<-100
phis<-seq(-5,10,length.out = C)
M<-10000
set.seed(99)
alpha_rate<-rep(0,C)
for (j in 1:C) {
  counter <-0
  for (i in 2:M) {
    y<-rnorm(1,xis[i-1],sqrt(exp(phis[j])))
    alph<-min(1,exp(theta*y-exp(theta*y)-(theta*xis[i-1]-exp(theta*xis[i-1]))))
    u<-runif(1)
    if (alph >= u){
      counter<-counter+1
      xis[i]<-y
    }
    else{
      xis[i]<-xis[i-1]
    }
  }
  alpha_rate[j]<-counter/M
}
plot(phiC,acceptRate)
axis(2,c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
axis(1,c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10))

j=1
while(alpha_rate[j]>.234){
  j=j+1
}

if(abs(alpha_rate[j-1]-.234)>abs(alpha_rate[j]-.234)){
  phis[j]} else {phis[j-1]}