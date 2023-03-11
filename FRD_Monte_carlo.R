set.seed(1221126)
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
library(MuMIn)

# x<-runif(100,-2,2)
# ran<- sample(0:100, 10, TRUE)
# y<-ifelse(x>=0,x+10,x+1)#from here we know that the treatment affect should be 9
# 
# for (i in ran){
#   if (x[i] <= 0){
#     y[i] = 3*x[i]+10
#   }
#   else {y[i] = 3*x[i]+1}
# }
# plot(x, y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
# abline(v = 0)

#Data generating process With always taker and never taker
#Linear_FRD
DGF1<-function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x+1+gap,x+1)
  y<-y+rnorm(n,0,1)
  w_s<-ifelse(x>=0,1,0)
  ran<- sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]+1+rnorm(1,0,1),
                   x[i]+1+gap+rnorm(1,0,1))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}

#Non_Linear_FRD
DGF2 <- function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+1+gap,-x^2+1)
  y<-y+rnorm(n,0,5)
  w_s<-ifelse(x>=0,1,0)
  ran<-sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+1+rnorm(1,0,5),
                   -x[i]^2+1+gap+rnorm(1,0,5))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}

#Third_order_FRD
DGF3 <- function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+gap,x^2+x^3+1)
  y<-y+rnorm(n,0,20)
  w_s<-ifelse(x>=0,1,0)
  ran<- sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+x[i]^3+1+rnorm(1,0,20),
                   x[i]^2+x[i]^3+gap+rnorm(1,0,20))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}

#Fourth_order_FRD
DGF4<-function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+x^4+1+gap,-x^2+x^3-x^4+1)
  y<-y+rnorm(n,0,50)
  w_s<-ifelse(x>=0,1,0)
  ran<- sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+x[i]^3+x[i]^4+1+rnorm(1,0,50),
                   -x[i]^2+x[i]^3-x[i]^4+1+gap+rnorm(1,0,50))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}


# #Checking the data generating process
data <- DGF1(100,10,9)
plot(data$x, data$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
     main="DGF1(n=100,a=10,gap=9)")
abline(v = 0)
# loc<-rdrobust(data$y,data$x,fuzzy=data$w_s,p=1,c=0,kernel="uniform")
# loc$bws[1,]
data2 <- DGF2(100,10,99)
plot(data2$x, data2$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
    main="DGF2(n=100,a=10,gap=99)" )
abline(v = 0)
# 
# data3 <- DGF3(100,10,199)
# plot(data3$x, data3$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
#      main="DGF3(n=100,a=10,gap=199)")
# abline(v = 0)
# 
# data4 <- DGF4(100,10,999)
# plot(data4$x, data4$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
#      main="DGF4(n=100,a=10,gap=999)")
# abline(v = 0)

#Preparation for Monte Carlo
kernel <- c("triangular","epanechnikov","uniform")

#Monte Carlo simulation
FRD <- function(T,order,n,DGP,kernel,gap){
  mse_kernel <- matrix(NA,length(kernel),order)  
  for (ker in 1:length(kernel)){
    a<-c()
    np<-matrix(NA,T,order)
    for (i in 1:T){
      for (k in 1:order){
        a[i]<-20
        data<-DGP(n,a[i],gap)
        locfit<-rdrobust(data$y,data$x,fuzzy=data$w_s,p=k,c=0,kernel=kernel[ker])
        np[i,k]<-locfit$coef[1]
      }
  }
    for (z in 1:order){
    mse_kernel[ker,z] <- mean((np[1:T,z]-gap)^2)}
      #sqrt(mean((np[1:T,z]-gap)^2))/mean(np[1:T,z])}
  }
return(mse_kernel)
}

#DGP can be DGF1,DGF2,DGF3 and DGF4, which is for non-linear simulation
mse_f1<- FRD(1000,4,DGP=DGF1,kernel=kernel,n=1000,gap = 9)#n cannot be too small
mse_f2<- FRD(1000,4,DGP=DGF2,kernel=kernel,n=1000,gap = 9)
mse_f3<- FRD(1000,4,DGP=DGF3,kernel=kernel,n=1000,gap = 9)
mse_f4<- FRD(1000,4,DGP=DGF4,kernel=kernel,n=1000,gap = 9)


###########different polynomial :
n_f<-seq(1000,10000,by=1000)
p=4
# mse_nf_t<-matrix(NA,length(n_f),4)
mse_nf_e<-matrix(NA,length(n_f),4)
# mse_nf_u<-matrix(NA,length(n_f),4)
for (i in 1:length(n_f)){
  # mse_nf_t[i,]<-FRD(1000,p,n_f[i],DGF1,kernel[1],gap=9)
  mse_nf_e[i,]<-FRD(1000,p,n_f[i],DGF1,kernel[2],gap=9)
  # mse_nf_u[i,]<-FRD(1000,p,n_f[i],DGF1,kernel[3],gap=9)
}
mse_nf_e
#then change the argument of mse_nf_u: DGF1 to DGF2,DGF3,DGF4.
#so Then we can get all the MSE of different n when kernel is "epanechnikov".






########find the critical value of the gap:
gap_s<-c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001,0.000000001,0.0000000001)
##for DGF1:
m_f2_1<-matrix(NA,length(gap_s),4)
m_f2_2<-matrix(NA,length(gap_s),4)
m_f2_3<-matrix(NA,length(gap_s),4)
m_f2_4<-matrix(NA,length(gap_s),4)
for (i in 1:length(gap_s)){
  m_f2_4[i,]=FRD(1000,4,DGP=DGF4,kernel=kernel[2],n=1000,gap = gap_s[i])#n cannot be too small
}


###Parametric Method
vars <- c("w_s","data$x_centered","I(x_centered^2)","I(x_centered^3)","I(x_centered^4)")
para <- function(T,vars,c,DGP,gap,n){
  a <- c()
  coe <- matrix(NA,T,4)
  for (i in 1:T){
    a[i]<-runif(1,min=1,max=10)
    data<-DGP(n,a[i],gap)
    data$x_centered<-data$x-c
    for (z in 1:(length(vars)-1)){
      formula <- as.formula(paste("y",paste(vars[0:z+1],collapse = "+"),sep = "~"))
      fit <- lm(formula,data=data)
      coe[i,z] <- fit$coefficients[2]
    }
  }
  sqr <- (coe-gap)^2
  mse <- apply(sqr,2,mean)
  return(mse)
}
##DGP can be DGF or DGF2 or DGF3, which is for non-linear data generating process
avg_f<-matrix(NA,4,4) 
avg_f[1,]<- para(T,vars=vars,c=0,DGP=DGF1,gap=9,n=1000)
avg_f[2,]<- para(T,vars=vars,c=0,DGP=DGF2,gap=9,n=1000)
avg_f[3,]<- para(T,vars=vars,c=0,DGP=DGF3,gap=9,n=1000)
avg_f[4,]<- para(T,vars=vars,c=0,DGP=DGF4,gap=9,n=1000)

#Simple testing
data <- DGF3(100,30)
fit <- lm(y~w_s+x,data=data)
fit$coefficients

ma <- matrix(4,2,2)
q <- ma-1
q_2 <- q^2
q_2

