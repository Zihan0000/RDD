rm(list=ls())
set.seed(122122)
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
library(MuMIn)
DGF1<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x+10,x+1)
  y<-y+rnorm(n,0,1)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
DGF2 <- function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+10,x^2+1)
  y<-y+rnorm(n,0,5)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
DGF3<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+10,x^2+x^3+1)
  y<-y+rnorm(n,0,20)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
DGF4<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+x^4+10,x^2+x^3+x^4+1)
  y<-y+rnorm(n,0,50)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
##non-parametric method:
kernel <- c("triangular","epanechnikov","uniform")
## we will see the treatment effect in different kernel:
FRD <- function(T,order,n,DGP,kernel){
  mean_kernel <- matrix(NA,length(kernel),order)  
  for (ker in 1:length(kernel)){
    a<-c()
    np<-matrix(NA,T,order)
    for (i in 1:T){
      for (k in 1:order){
        a[i]<-runif(1,min=1,max=20)
        data<-DGP(n,a[i])
        locfit<-rdrobust(data$y,data$x,p=k,c=0,kernel=kernel[ker])
        np[i,k]<-locfit$coef[1]
      }
    }
    for (z in 1:order){
      mean_kernel[ker,z] <- mean(np[1:T,z])}
  }
  return(mean_kernel)
}
#then we can use the above fuction to calculate the treatment effect in different kernel:
T<-1000
p<-4
n<-300
#I hope I can think a way to make the below codes simpler.
mean1<-matrix(NA,3,p)#mean1 is the mean value of the treatment effect in 1-polynomial 
for (i in 1:3){
  mean1[i,]<-FRD(T,p,n,DGF1,kernel[i])
}
mean2<-matrix(NA,3,p)#mean2 is the mean value of the treatment effect in 2-polynomial 
for (i in 1:3){
  mean2[i,]<-FRD(T,p,n,DGF2,kernel[i])
}
mean3<-matrix(NA,3,p)#mean3 is the mean value of the treatment effect in 3-polynomial 
for (i in 1:3){
  mean3[i,]<-FRD(T,p,n,DGF3,kernel[i])
}
mean4<-matrix(NA,3,p)#mean4 is the mean value of the treatment effect in 1-polynomial 
for (i in 1:3){
  mean4[i,]<-FRD(T,p,n,DGF4,kernel[i])
}
 #parametric method:
vars <- c("w_s","data$x_centered","I(data$x_centered^2)","I(data$x_centered^3)","I(data$x_centered^4)")
para <- function(T,vars,c,DGP){
  a<-c()
  coe <- matrix(NA,T,4)
  for (i in 1:T){
    a[i]<-runif(1,min=1,max=20)
    data<-DGP(n,a[i])
    data$x_centered<-data$x-c
    for (z in 1:(length(vars)-1)){
      formula <- as.formula(paste("y",paste(vars[0:z+1],collapse = "+"),sep = "~"))
      fit <- lm(formula,data=data)
      coe[i,z] <- fit$coefficients[2]
    }
  }
  avg <- apply(coe, 2, mean)
  return(avg)
}
avg<-matrix(NA,4,4)
avg[1,]<- para(T,vars=vars,c=0,DGP=DGF1)
avg[2,]<- para(T,vars=vars,c=0,DGP=DGF2)
avg[3,]<- para(T,vars=vars,c=0,DGP=DGF3)
avg[4,]<- para(T,vars=vars,c=0,DGP=DGF4)


