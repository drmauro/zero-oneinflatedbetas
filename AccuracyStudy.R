
rm(list=ls(all=TRUE))
#set.seed(0)

start_time <- Sys.time()


library(magrittr)

library(foreach)

library(doParallel)


#setwd('C:\\Users\\jafio\\Dropbox\\Mistura Betas Inflacionadas')
setwd('C:\\Users\\Augusto\\Dropbox\\Mistura Betas Inflacionadas')


MyData <- read.csv(file="dados.csv", header=TRUE, sep=",")
summary(MyData)


X <- MyData[, !(names(MyData) %in% "lgd")] %>% scale
k = dim(X)[2]
X2 = data.frame(inter = rep(1, length(X[,1])) , X)
X_m = X2 %>% as.matrix

y <- MyData[,  (names(MyData) %in% "lgd")]



################################################################
#### Inflated Mixed Beta Model with 0 and with 1   #############
################################################################

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


### SSE: sum of squares of errors
##  regressão para todos os parâmetros
SSE_reg <- function(varp, y, X_m){
  
  lgd_hat <- y.fitted(varp, X_m) 
  
  sum( (y - lgd_hat)^2 )
}

################################################################
################################################################



################################################################
#### Modelo Beta Inflacionado no 0 e no 1       ################
################################################################
y.fitted_reg_beta01 <- function(varp, X_m){
  varp_m = matrix(varp, ncol = 4, byrow = F)
  matrix_pl = X_m %*% varp_m
  
  pl1 = matrix_pl[,1]
  pl2 = matrix_pl[,2]
  #pl3 = matrix_pl[,3]
  pl4 = matrix_pl[,3]
  pl5 = matrix_pl[,4]
  #pl6 = matrix_pl[,6]
  #pl7 = matrix_pl[,7]
  
  
  gamm_zer = exp(pl1)/(1+exp(pl1)+exp(pl2))  
  gamm_one = exp(pl2)/(1+exp(pl1)+exp(pl2))
  
  #Pii      =  exp(pl3)/(1+exp(pl3))
  
  mu1      = exp(pl4)/(1+exp(pl4))  
  phi1     = exp(pl5)
  
  #mu2      = exp(pl6)/(1+exp(pl6))  
  #phi2     = exp(pl7)
  
  #lgd_hat = (1-gamm_zer-gamm_one)*(Pii*mu1+(1-Pii)*mu2) + gamm_one
  lgd_hat = (1-gamm_zer-gamm_one)*mu1 + gamm_one
  
  return(lgd_hat)
}

SSE_reg_beta01 <- function(varp, y, X_m){
  lgd_hat <- y.fitted_reg_beta01(varp, X_m) 
  
  sum( (y - lgd_hat)^2 )
}

SSE_reg_beta01(rep(0,36),y,X_m)

################################################################
################################################################




############## Estudo da capacidade preditiva do modelo #######
######## e comparação com a Beta Inflacionada #################
NSIM <- 1000
prop <- 0.7

ini <- c(0.572,0.227,-0.127,1.714,-1.522,-0.088,0.068,
         0.183,-1.268,0.083,0.083,0.429,-5.43,0.599,
         0.564,-0.042,-0.249,0.874,0.217,0.041,-0.102,
         -0.582,-0.006,0.037,-0.144,0.308,-1.242,1.245,
         0.161,-0.584,-0.344,1.868,-0.101,2.775,-0.192,
         0.328,0.109,0.421,-0.413,-0.087,1.009,-0.412,
         -1.044,-0.987,-0.059,1.308,0.492,-0.216,-0.745,
         0.433,-2.362,-0.064,-0.922,0.225,-0.467,1.589,
         -1.587,1.436,0.118,0.18,0.263,1.164,-1.042)


n = length(MyData[,1])

nn <- floor(prop*n)

matrix_par     <- matrix(NA, nrow = NSIM , ncol = 63 )
matrix_par.B     <- matrix(NA, nrow = NSIM , ncol = 36 )
matrix_idx.fit <- matrix(NA, nrow = NSIM , ncol = n )
matrix_y.test <-  matrix(NA, nrow = NSIM , ncol = n-nn )
matrix_y.hat  <-  matrix(NA, nrow = NSIM , ncol = n-nn )
matrix_y.hat.B  <-  matrix(NA, nrow = NSIM , ncol = n-nn )


count = 0
while (count < NSIM) {
  
  count = count+1
  print(count)
  
  idx.fit <- 1:n %in% sample.int(n=n, size=nn, replace=F)
  
  X.fit = X_m[idx.fit, ]
  y.fit = y[idx.fit ]
  
  X.test = X_m[!idx.fit, ]
  y.test = y[!idx.fit ]
  
  matrix_idx.fit[count, ] <- idx.fit
  
  matrix_y.test[count, ] <- y.test
  
}




registerDoParallel( detectCores()-1 )       ### Inicio da PARTE que RODA em PARALELO

####### Modelo A: Mistura de Betas Inflacionado no 0 e no 1 #######################
matrix_par = foreach (count=1:NSIM, .combine='rbind') %dopar% {  ## 'for' em paralelo 
  
  idx.fit <- matrix_idx.fit[count, ] 
  
  X.fit = X_m[idx.fit, ]
  y.fit = y[idx.fit ]
  
  EMQ = optim(ini, SSE_reg, method="Nelder-Mead", control = list(maxit=10^4, reltol=0.01), hessian=F, y=y.fit, X_m=X.fit )
  
  EMQ$par
}

matrix_y.hat = foreach (count=1:NSIM, .combine='rbind') %dopar% {
  idx.fit <- matrix_idx.fit[count, ]
  X.test  <- X_m[!idx.fit, ]
  varp    <- matrix_par[count, ]
  y.fitted(varp, X_m=X.test)
}
###################################################################################


####### Modelo B: Beta Inflacionado no 0 e no 1 ###################################
matrix_par.B = foreach (count=1:NSIM, .combine='rbind') %dopar% {  ## 'for' em paralelo 
  
  idx.fit <- matrix_idx.fit[count, ] 
  
  X.fit = X_m[idx.fit, ]
  y.fit = y[idx.fit ]
  
  EMQ_B = optim(rep(0,36), SSE_reg_beta01, method="Nelder-Mead", control = list(maxit=10^4, reltol=0.01), hessian=F, y=y.fit, X_m=X.fit )
  
  EMQ_B$par
}

matrix_y.hat.B = foreach (count=1:NSIM, .combine='rbind') %dopar% {
  idx.fit <- matrix_idx.fit[count, ]
  X.test  <- X_m[!idx.fit, ]
  varp    <- matrix_par.B[count, ]
  y.fitted_reg_beta01(varp, X_m=X.test)
}
###################################################################################


stopImplicitCluster()                         ### Fim da PARTE que RODA em PARALELO







######## Métricas ################################################################

tab <- matrix(NA, 8, 2)
colnames(tab) <- c('Mix Beta', 'Beta')
rownames(tab) <- c('Mean Error', 'Std. Desv. Error', 'Q1 Error ', 'Q2 Error', 'Q3 Error', 'RMSE', 'MSE', 'MAE')

### Modelo A: modelo de mistura de 2 Betas Inflacionado

eA <- na.omit(matrix_y.test - matrix_y.hat) # Erro

tab[1, 1] = mean(eA)
tab[2, 1] = sd(eA)
tab[3:5, 1] = quantile(eA, c(0.25, 0.5, 0.75))

# RMSE
tab[6, 1] = sqrt( mean( (matrix_y.test - matrix_y.hat)^2 , na.rm=T ) )
# MSE
tab[7, 1] = mean( (matrix_y.test - matrix_y.hat)^2, na.rm=T  ) 
# MAE
tab[8, 1] = mean( abs(matrix_y.test - matrix_y.hat), na.rm=T  ) 




### Modelo B: modelo Beta Inflacionado

eB <- na.omit(matrix_y.test - matrix_y.hat.B) # Erro

tab[1, 2] = mean(eB)
tab[2, 2] = sd(eB)
tab[3:5, 2] = quantile(eB, c(0.25, 0.5, 0.75))

# RMSE
tab[6, 2] = sqrt( mean( (matrix_y.test - matrix_y.hat.B)^2, na.rm=T ) )
# MSE
tab[7, 2] = mean( (matrix_y.test - matrix_y.hat.B)^2, na.rm=T  ) 
# MAE
tab[8, 2] = mean( abs(matrix_y.test - matrix_y.hat.B), na.rm=T  ) 

tab <- tab %>% round(4)

t(tab)

##  write.csv(tab, 'Estudo_Acuracia.csv')

###################################################################################


end_time <- Sys.time()

## tempo total
end_time - start_time

## save.image("EstudoAcuracia_1000.Rdata")


###################################################################################
########### Histogramas dos parametros #############################################
###################################################################################

load("EstudoAcuracia_1000.Rdata")

## gamm_zer
par(mfrow=c(3,3))

for(i in 1:9)
  hist( matrix_par[,i], main= paste0('beta_1,',i) , xlab='')


## gamm_one
par(mfrow=c(3,3))

for(i in 10:18)
  hist( matrix_par[,i], main= paste0('beta_2,',i) , xlab='')


## Pii
par(mfrow=c(3,3))

for(i in 19:27)
  hist( matrix_par[,i], main= paste0('beta_3,',i) , xlab='')


## mu1
par(mfrow=c(3,3))

for(i in 28:36)
  hist( matrix_par[,i], main= paste0('beta_4,',i) , xlab='')


## phi1
par(mfrow=c(3,3))

for(i in 37:45)
  hist( matrix_par[,i], main= paste0('beta_5,',i) , xlab='')


## mu2
par(mfrow=c(3,3))

for(i in 46:54)
  hist( matrix_par[,i], main= paste0('beta_6,',i) , xlab='')


## phi2
par(mfrow=c(3,3))

for(i in 55:63)
  hist( matrix_par[,i], main= paste0('beta_7,',i) , xlab='')

###################################################################################

