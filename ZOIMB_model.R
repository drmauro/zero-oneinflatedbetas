rm(list=ls(all=TRUE))
#set.seed(0)

library(magrittr)
library(GA)

# Set the working diretory 
setwd('C:\\Users\\jafio\\Dropbox\\Mistura Betas Inflacionadas')
#setwd("C:/Users/Mauro/Dropbox/Mistura Betas Inflacionadas (1)")

MyData <- read.csv(file="dados.csv", header=TRUE, sep=",")
summary(MyData)

MyData %>% dim
MyData %>% head
MyData %>% names

#Install if the package doesn't exist 
#install.packages('DataExplorer) 
library(DataExplorer)
plot_str(MyData)
plot_missing(MyData)
plot_histogram(MyData)
plot_correlation(MyData, type = 'continuous','Review.Date')


## Ajuste Pontual
X = MyData[, !(names(MyData) %in% "lgd")]
X = scale(X)
y = MyData[,  (names(MyData) %in% "lgd")]

k = dim(X)[2]
X2 = data.frame(inter = rep(1, length(X[,1])) , X)
X_m = X2 %>% as.matrix
dim(X_m)
X_m %>% head



###############################################################################
#############        Zero One Inflated Mix Beta Model               ###########
###############################################################################
# density
dZOIMM.call <- function(y,gamma0,gamma1,pi,mu1,phi1,mu2,phi2){
  
  if(y==0)
    return(gamma0)
  if(y==1)
    return(gamma1)
  
  fb1  <- dbeta(y,mu1*phi1,(1-mu1)*phi1) #pdf Beta 1
  fb2  <- dbeta(y,mu2*phi2,(1-mu2)*phi2) #pdf Beta 2
  fM   <- (pi*fb1)+(1-pi)*fb2            #pdf Mix Beta 1 e Beta 2
  return( (1-(gamma0+gamma1))*fM )
}

dZOIMM <- function(y,gamma0,gamma1,pi,mu1,phi1,mu2,phi2){

  if(gamma0+gamma1 >1 | gamma0 <0 | gamma1 <0 | pi <0 | pi >1 | mu1 <0 | mu1 >1 | mu2 <0 | mu2 >1 | phi1 <0 | phi2 <0 | mu1 > mu2) 
    return( rep(0,length(y)) )
  
  sapply(y, function(y) dZOIMM.call(y,gamma0,gamma1,pi,mu1,phi1,mu2,phi2) )
}


# cdf - cumulative distribution function
pZOIMM.call <- function(y,gamma0,gamma1,pi,mu1,phi1,mu2,phi2){
  if(y < 0)
    return(0)
  if(y == 0)
    return(gamma0)
  if(y >= 1)
    return(1)
  
  pb1  <- pbeta(y,mu1*phi1,(1-mu1)*phi1) #cdf Beta 1
  pb2  <- pbeta(y,mu2*phi2,(1-mu2)*phi2) #cdf Beta 2
  pM   <- (pi*pb1)+(1-pi)*pb2            #cdf Mix Beta 1 e Beta 2
  
  return( gamma0 + (1-gamma0-gamma1)*pM )
}

pZOIMM <- function(y,gamma0,gamma1,pi,mu1,phi1,mu2,phi2){
  sapply(y, function(y) pZOIMM.call(y,gamma0,gamma1,pi,mu1,phi1,mu2,phi2) )
}

#log likelihood function
loglike <- function(par, y){
  dZOIMM(y,par[1],par[2],par[3],par[4],par[5],par[6],par[7]) %>% log %>% sum
}
###############################################################################
###############################################################################


###############################################################################
############# Zero One Inflated Beta Model (Ospina & Ferrari, 2010) ###########
###############################################################################
dZOIM <- function(y,gamma0,gamma1,mu1,phi1){
  
  if(gamma0+gamma1 >1 | gamma0 <0 | gamma1 <0 | mu1 <0 | mu1 >1 | phi1 <0 ) 
    return( rep(0,length(y)) )
  
  sapply(y, function(y) dZOIMM.call(y,gamma0,gamma1,pi=1,mu1,phi1,mu2=0.99,phi2=1) )
}

pZOIM <- function(y,gamma0,gamma1,mu1,phi1){
  sapply(y, function(y) pZOIMM.call(y,gamma0,gamma1,pi=1,mu1,phi1,mu2=0.99,phi2=1) )
}

loglike_ZOIM <- function(par, y){
  dZOIMM(y,par[1],par[2],1,par[3],par[4],0.99,1) %>% log %>% sum
}
###############################################################################
###############################################################################



#############################################################################
###  DATA WITHOUT ZERO's AND WITHOUT ONE's
#############################################################################

z <- y[ !(y==0 | y==1) ]
hist(z,100)

#### Mix Beta Model
loglikeMM <- function(par, z){
  loglike( par=c(0,0,par), y=z)
}

#opt.z <- optim(c(0.5, 0.3, 1, 0.7, 1), loglikeMM, control= list(fnscale=-1, maxit=10^4),z=z)
#opt.z$par
#[1]  0.26770519  0.07649713 12.17261162  0.86036774  9.35425086

yy <- seq(0,1,0.01)
dd <- dZOIMM(yy, gamma0=0, gamma1=0, pi=0.26770519, mu1=0.07649713, phi1=12.17261162, mu2=0.86036774, phi2=9.35425086)

par(mfrow=c(1,1))
hist(z,100,freq =F)
points(yy,dd,type='l', col='red')


#############################################################################
########### Maximum Likelihood Estimation   #################################
#############################################################################


## ZOIMB ##############
par_ini <- c(mean(y==0), mean(y==1), 0.26770519,  0.07649713, 12.17261162,  0.86036774,  9.35425086)
opt <- optim(par_ini, loglike, hessian = TRUE, control= list(fnscale=-1, maxit=10^4),y=y)
opt$par
opt$hessian %>% solve %>% diag %>% abs %>% sqrt %>% round(4)

## Best optimizer option
# opt.ga <- ga(type="real-valued", fitness=loglike, lower=c(0,0,0,0,0,0,0), upper=c(0.5,0.5,1,0.9,10,0.9,10), y=y)


## AIC
-2*loglike(opt$par, y) + 2*7


yy <- seq(-0.2,1.2,0.01)
dd <- dZOIMM(yy, gamma0=0.32372177, gamma1=0.30901536, pi=0.26770519, mu1=0.07649713, phi1=12.17261162, mu2=0.86036774, phi2=9.35425086)
pp <- pZOIMM(yy, gamma0=0.32372177, gamma1=0.30901536, pi=0.26770519, mu1=0.07649713, phi1=12.17261162, mu2=0.86036774, phi2=9.35425086)


## ZOIB ou ZOIM ####################
opt_ZOIM <- optim(c(0.3, 0.3, 0.7, 1), loglike_ZOIM, hessian = TRUE, control= list(fnscale=-1, maxit=10^4),y=y)
opt_ZOIM$par

opt_ZOIM$hessian %>% solve %>% diag %>% abs %>% sqrt %>% round(4)


## AIC
-2*loglike_ZOIM(opt_ZOIM$par, y) + 2*4


yy <- seq(-0.2,1.2,0.01)
dd_ZOIM <- dZOIM(yy, gamma0=0.3237445, gamma1=0.3089484, mu1=0.5835252, phi1=1.1256703)
pp_ZOIM <- pZOIM(yy, gamma0=0.3237445, gamma1=0.3089484, mu1=0.5835252, phi1=1.1256703)





### parece não fazer sentido esse gráfico nesse caso
par(mfrow=c(1,1))
#hist(y, seq(0,1,0.1), probability=T)
hist(y, 100, probability=T)
points(yy,dd,type='l', col='red')


### Distribuição acumulada
par(mfrow=c(1,1))
ecdf(y) %>% plot(verticals = TRUE, do.points = FALSE, main='', xlab='y', ylab='F( y )', col='gray', lwd=5)
points(yy,pp,type='l', col='red', lwd=3, lty=2)
points(yy,pp_ZOIM,type='l', col='blue', lwd=3 , lty=3)

legend(x=0.7, y=0.3, legend=c('Empirical', 'ZOIMB','ZOIB'), col=c(8,2,4), lwd=c(5,2,2),lty=c(1,2,3))


########################################################################
#########   ZOIMM with regressor variables     #########################
########################################################################

y.fitted <- function(varp, X_m){
  varp_m = matrix(varp, ncol = 7, byrow = F)
  matrix_pl = X_m %*% varp_m
  
  pl1 = matrix_pl[,1]
  pl2 = matrix_pl[,2]
  pl3 = matrix_pl[,3]
  pl4 = matrix_pl[,4]
  pl5 = matrix_pl[,5]
  pl6 = matrix_pl[,6]
  pl7 = matrix_pl[,7]
  
  
  gamm_zer = exp(pl1)/(1+exp(pl1)+exp(pl2))  
  gamm_one = exp(pl2)/(1+exp(pl1)+exp(pl2))
  
  Pii      =  exp(pl3)/(1+exp(pl3))
  
  mu1      = exp(pl4)/(1+exp(pl4))  
  phi1     = exp(pl5)
  
  mu2      = exp(pl6)/(1+exp(pl6))  
  phi2     = exp(pl7)
  
  lgd_hat = (1-gamm_zer-gamm_one)*(Pii*mu1+(1-Pii)*mu2) + gamm_one
  
  return(lgd_hat)
}


### soma dos quadrados dos erros: sum of squares of errors
##  regressão para todos os parâmetros: regression for all parameters
SQE_reg <- function(varp, y, X_m){
  
  lgd_hat <- y.fitted(varp, X_m) 
  
  sum( (y - lgd_hat)^2 )
}

#ini_0 <- rep(0,63)
#EMQ = optim(ini_0, SQE_reg, method="Nelder-Mead", control = list(maxit=10^4, reltol=0.01), hessian=F, y=y, X_m=X_m )

#ini <- EMQ$par
ini <- c(0.572,0.227,-0.127,1.714,-1.522,-0.088,0.068,
         0.183,-1.268,0.083,0.083,0.429,-5.43,0.599,
         0.564,-0.042,-0.249,0.874,0.217,0.041,-0.102,
         -0.582,-0.006,0.037,-0.144,0.308,-1.242,1.245,
         0.161,-0.584,-0.344,1.868,-0.101,2.775,-0.192,
         0.328,0.109,0.421,-0.413,-0.087,1.009,-0.412,
         -1.044,-0.987,-0.059,1.308,0.492,-0.216,-0.745,
         0.433,-2.362,-0.064,-0.922,0.225,-0.467,1.589,
         -1.587,1.436,0.118,0.18,0.263,1.164,-1.042)

EMQ = optim(ini, SQE_reg, method="Nelder-Mead", control = list(maxit=10^4, reltol=0.01), 
            hessian=F, y=y, X_m=X_m )

y.hat <- y.fitted( EMQ$par, X_m)

hist(y.hat,100, main='')

residuo <- y - y.hat

residuo %>% summary

plot(residuo)

hist(residuo,50, main='')

plot(y, y.hat)


########################################################################
########################################################################




