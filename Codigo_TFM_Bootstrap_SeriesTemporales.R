#########################################################################
## Código para obtener los resultados del TFM 
## Bootstrap para series temporales en modelos de producción excedentaria
#########################################################################

### Librerias a utilizar
library(mgcv)
library(fpp2) 
library(tseries)
library(ggplot2)
### Paquete para utilizar knobi
### El primer paso es instalar el paquete devtools usando 
# install.packages("devtools")

### Hay dos opciones para instalar el paquete:
### Opción 1, versión con vignettes
# install.packages(c("corrplot", "ggplot2",  "gridExtra", "grDevices",  "optimr", "plot3D"))
# devtools::install_github("MERVEX-group/knobi",build_vignettes = TRUE)

### Opción 2, versión sin vignettes
devtools::install_github("MERVEX-group/knobi")
### Esta es una opción más rápida.
### Se recomienda reiniciar R después de la instalación del paquete.
### libreria
library(knobi)
### En caso de tener problemas con knobi debido a la versión del ggplot2
### knobi actualmente funciona con la versión siguiente
# devtools::install_version("gtable", version = "0.3.1", repos = "http://cran.us.r-project.org")

### Inicio del código

### Semilla utilizada
set.seed(123)
### Número de simulaciones: réplicas Montecarlo
nsim=200
### Tamaño de la serie
tmax=100

### Parámetros reales
real_p=1.2
real_K=110.8215
real_r=0.46441
real_Bmsy <- real_K*(1/(real_p+1))^(1/real_p)
real_MSY <- real_r*real_K*(1/(real_p+1))^(1+1/real_p)
real_Fmsy = (real_r/real_p)*(1-(1/(real_p+1)))

### Función de simulación Montecarlo
Simulation=function(p, r, K, nsim=100, tmax=100){
  
  p=p
  K=K
  r=r
  
  tmax=tmax
  nsim=nsim
  
  dat_iter_C=matrix(0, ncol = tmax, nrow = nsim)
  dat_iter_B=matrix(0, ncol = tmax+1, nrow = nsim)
  dat_iter_F=dat_iter_C
  
  B0=K
  
  Bmsy <- K*(1/(p+1))^(1/p)
  MSY <- r*K*(1/(p+1))^(1+1/p)
  x <- MSY/K
  Fmsy = (r/p)*(1-(1/(p+1)))
  Fcrash = r/p
  
  mean_1 <- 0.5*Fmsy
  mean_2 <- 0.9*Fcrash
  x_approx <- c(mean_1,mean_2)
  data_approx1 <- approx(x_approx, n=tmax/5) 
  x_approx2 <- c(mean_2,Fmsy)
  data_approx2 <- approx(x_approx2, n=tmax/5)
  x_approx3 <- c(0.25*Fmsy,mean_1)
  data_approx3 <- approx(x_approx3, n=tmax/10) 
  
  Fmean <- c(0,data_approx3$y, rep(mean_1, tmax/5-1), data_approx1$y, rep(0.9*Fcrash,tmax/5), data_approx2$y, rep(Fmsy, tmax/10))
  
  cv <- 0.1
  Fsd <- cv * Fmean
  Ftrends <- mapply(rnorm,nsim,Fmean,Fsd)
  Ftrends[Ftrends<0] <- 0
  Ftrends

  for (i in 1:nsim) {
    
    C <- numeric()
    B <- numeric()
    B[1] = B0
    
    
    Fvec=Ftrends[i,] # Vector de diferentes mortalidades por pesca
    
    for(t in 1:tmax) {
      C[t] <- Fvec[t]*B[t]
      B[t+1] <- B[t]+(r/p)*B[t]*(1-(B[t]/K)^p)-C[t] 
      
    }
    
    
    dat_iter_C[i,] <- C
    dat_iter_B[i,] <- B 
    dat_iter_F[i,] <- Fvec 
  }
  
  dat_iter=list(catch=dat_iter_C,
                ssb=dat_iter_B,
                f=dat_iter_F)
  
  return ( dat_iter)
}

### Simulación de nuestro procedimiento Montecarlo 
simu<-Simulation(p=real_p, r=real_r, K=real_K, nsim=nsim, tmax=tmax);#head(simuB)
### Datos de biomasa
simuB<-simu[["ssb"]]

### Gráfico de las réplicas Montecarlo de biomasa
time<-1:tmax
btime=1:(tmax+1)

q <- matrix(0, ncol= 3, nrow= tmax+1) 

for(t in 1:(tmax+1)) {
  q[t,] <- quantile(as.numeric(simuB[,t]), c(0.025, 0.5, 0.975))
}

plot(x=1:(tmax+1),y=simuB[1,],type="l",col="gray",
     main = "Réplicas Montecarlo de biomasa",
     xlab="Time", ylab="Réplicas",
     ylim=c(min(simuB)*0.5,max(simuB)*1.1))

for (i in 2:nsim) {
  lines(btime, simuB[i,],type="l",col="gray")
}

lines(btime, q[,1], col="blue",lty=2)
lines(btime, q[,2], col="blue",lty=2)
lines(btime, q[,3], col="blue",lty=2)
  

### Datos de capturas 
simuC<-simu[["catch"]]

### Gráfico de las réplicas Montecarlo de capturas
q <- matrix(0, ncol= 3, nrow= tmax)

for(t in 1:tmax) {
  q[t,] <- quantile(as.numeric(simuC[,t]), c(0.025, 0.5, 0.975))
}

plot(time,simuC[1,],type="l",col="gray",
     main = "Réplicas Montecarlo de capturas",
     xlab="Time", ylab="Réplicas",
     ylim=c(min(simuC)*0.5,max(simuC)*1.1))

for (i in 2:nsim) {
  lines(time, simuC[i,],type="l",col="gray")
}

lines(time, q[,1], col="blue",lty=2)
lines(time, q[,2], col="blue",lty=2)
lines(time, q[,3], col="blue",lty=2)

### Datos de la mortalidad por pesca (F; Fishing)
simuF<-simu[["f"]]

### Gráfico de las réplicas Montecarlo de Fishing 
q <- matrix(0, ncol= 3, nrow= tmax)

for(t in 1:tmax) {
  q[t,] <- quantile(as.numeric(simuF[,t]), c(0.025, 0.5, 0.975))
}

plot(time,simuF[1,],type="l",col="gray",
     main = "Réplicas Montecarlo de fishing (F)",
     xlab="Time", ylab="Réplicas",
     ylim=c(min(simuF)*0.5,max(simuF)*1.1))

for (i in 2:nsim) {
  lines(time, simuF[i,],type="l",col="gray")
}

lines(time, q[,1], col="blue",lty=2)
lines(time, q[,2], col="blue",lty=2)
lines(time, q[,3], col="blue",lty=2)


### Obtenemos la biomasa a principio de año

simuB_average=simuB

for(i in 1:nsim){
  serie=simuB_average[i,]
  regresion=smooth.spline(serie~btime)
  simuB[i,]=predict(regresion,x=btime-0.5)$y
}

## Eliminamos un año al principio de cada serie
data <- list(catch=simuC[,-1], ssb=simuB[,-1])


## Función para knobi
knobi_fit_function=function(data, nsim=100){
  
  dat_par=matrix(0, ncol = 6, nrow = nsim);dat_par=data.frame(dat_par)
  colnames(dat_par) <- c("r", "p", "K", "B_MSY", "F_MSY", "MSY" )

  r <- numeric()
  p <- numeric()
  K <- numeric()
  B_MSY <- numeric() 
  F_MSY <- numeric()
  MSY <- numeric()
  
  for (i in 1:nsim) {
    
    knobi_data <- list( "Catch" = data[["catch"]][i,] , 
                  "Spawning_Biomass"= data[["ssb"]][i,])
    
    control<-list(pella="TRUE")
    
    # Fit the model
    knobi_results <-knobi::knobi_fit(data=knobi_data, control,plot_out = FALSE)
    
    if (any(seq(0,100,by=5)==100*(i/nsim))) {
      cat(paste0(100*(i/nsim)," % ... ","\n"))}
    
    r[i] <- knobi_results$fit$Parameter_estimates[1]
    p[i] <- knobi_results$fit$Parameter_estimates[3]
    K[i] <- knobi_results$fit$Parameter_estimates[2]
    B_MSY[i] <- knobi_results$fit$RP$B_MSY
    F_MSY[i] <- knobi_results$fit$RP$F_MSY
    MSY[i]  <- knobi_results$fit$RP$MSY
    
    dat_par[i,1] <- r[i]
    dat_par[i,2] <- p[i]
    dat_par[i,3] <- K[i]
    dat_par[i,4] <- B_MSY[i]
    dat_par[i,5] <- F_MSY[i]
    dat_par[i,6] <- MSY[i]
  
  }
  
  return ( dat_par)
}


### Estimaciónes de knobi (KBPMs) 
sim<-knobi_fit_function(data=data, nsim=nsim)
head(sim)


### Medidas de tendencia central y dispersión
### Media 
colMeans(sim)
summary(sim)

### Desviación estándar 
sapply(sim, sd)

### Cuantiles 5%, mediana y 95%

quantile_r <- quantile(sim$r, c(0.025, 0.5, 0.975))
quantile_p <- quantile(sim$p, c(0.025, 0.5, 0.975))
quantile_K <- quantile(sim$K, c(0.025, 0.5, 0.975))
quantile_B_MSY <- quantile(sim$B_MSY, c(0.025, 0.5, 0.975))
quantile_F_MSY <- quantile(sim$F_MSY, c(0.025, 0.5, 0.975))
quantile_MSY <- quantile(sim$MSY, c(0.025,  0.5, 0.975))

mc_comparison = rbind(quantile_r, quantile_p, quantile_K,
                              quantile_B_MSY, quantile_F_MSY, quantile_MSY)
rownames(mc_comparison) <- c("r", "p", "K", "B_MSY", "F_MSY", "MSY" )
mc_comparison <- cbind(mc_comparison,c(real_r,real_p,real_K,real_Bmsy,real_Fmsy,real_MSY))
colnames(mc_comparison)[4]="Real_values"
### Tabla comparativa de los intervalos de confianza y los parámetros reales
### de la simulación Montecarlo
mc_comparison


##########################################################################
### Aplicación al bootstrap
##########################################################################
### Número de remuestras por cada réplica Montecarlo
nboot=200

### Función para la obtención del bootstrap con GAM+AR
AR.gam.boot=function(data,nboot=100){
  
  vector=data
  
  time<-c(1:length(vector))
  n <- length(vector)
  gam_mod <- gam(vector~s(time))
  
  
  arma_res <- auto.arima(gam_mod$residuals, max.order=  1/4*(n),
                         max.p= 1/4*(n), max.q = 0, d=0,
                         stationary = TRUE,trace=FALSE,
                         seasonal = FALSE,
                         approximation = FALSE,
                         stepwise = TRUE, ic= "aic") 
  
  p<-arma_res$arma[1]
  
  if(p == 0){
    
    residuals_vector <- gam_mod$residuals
  
    } else {
    
    ar_mod<-Arima(gam_mod$residuals,order=c(p,0,0),include.mean = FALSE)
    
    residuals_vector <- ar_mod$residuals
  }
  
  
  dat_boot=matrix(0,ncol=n,nrow=(nboot))
  
  for(i in 1:nboot) {
    
    ind=sample(1:n,n, replace= FALSE)
    
    if(p == 0){
      
      fitted_values <- gam_mod$fitted.values 
    
      } else {
        
      fitted_values <- ar_mod$fitted + gam_mod$fitted.values 
    
    }
    
    dat_boot[i,] <- fitted_values + residuals_vector[ind]
    
  }
  
  q <- matrix(0, ncol= 3, nrow= n) 
  
  for(i in 1:n) {
    q[i,] <- quantile(dat_boot[,i], c(0.025, 0.5, 0.975))
  }
  
  # Gráfico de remuestras 
  plot(dat_boot[1,],type="l",col="gray",
       main = "Remuestras  con ajuste basado en modelos GAM",
       xlab="Time", ylab="Remuestras",
       ylim=c(min(dat_boot)*0.5,max(dat_boot)*1.1))
  
  for (i in 2:nboot) {
    lines(dat_boot[i,],type="l",col="gray")
  }
  
  lines(vector, type = "l",lwd=1 ,col= "black")
  lines(time, q[,1], col="blue",lty=2)
  lines(time, q[,2], col="blue",lty=2)
  lines(time, q[,3], col="blue",lty=2)
  
  return(dat_boot)
}


## Data para función knobi_fit_function del procedimiento Bootstrap

boot_c=boot_f=matrix(0,nrow=nboot*nsim,ncol=tmax-1)
boot_ssb=matrix(0,nrow=nboot*nsim,ncol=tmax)


for(i in 1:nsim){    
  
  for(j in 1:nboot){    

    boot_c[((i-1)*nboot+1):((i)*nboot),] = AR.gam.boot(data=data$catch[i,], nboot=nboot)
    boot_ssb[((i-1)*nboot+1):((i)*nboot),] = AR.gam.boot(data=data$ssb[i,], nboot=nboot)

  }
}

boot_data <- list(catch=boot_c, ssb=boot_ssb)


## Revisar las remuestras para ver si hay números negativos
which(boot_data$catch<0)
which(boot_data$ssb<0)

### Ajuste knobi_fit_function

simBoot<-knobi_fit_function(data=boot_data, nsim=nsim*nboot); head(simBoot)
head(simBoot)

### Medidas de tendencia central y dispersión
### Media
colMeans(simBoot)
summary(simBoot)

### Desviación estándar
sapply(simBoot, sd)

### Cuantiles 5%, mediana y 95%

quantile_r <- quantile(simBoot$r, c(0.025,0.5, 0.975))
quantile_p <- quantile(simBoot$p, c(0.025, 0.5, 0.975))
quantile_K <- quantile(simBoot$K, c(0.025, 0.5, 0.975))
quantile_B_MSY <- quantile(simBoot$B_MSY, c(0.025, 0.5, 0.975))
quantile_F_MSY <- quantile(simBoot$F_MSY, c(0.025, 0.5, 0.975))
quantile_MSY <- quantile(simBoot$MSY, c(0.025, 0.5, 0.975))

boot_comparison = rbind(quantile_r, quantile_p, quantile_K,
                      quantile_B_MSY, quantile_F_MSY, quantile_MSY)
rownames(boot_comparison) <- c("r", "p", "K", "B_MSY", "F_MSY", "MSY" )
boot_comparison <- cbind(boot_comparison,c(real_r,real_p,real_K,real_Bmsy,real_Fmsy,real_MSY))
colnames(boot_comparison)[4]="Real_values"

### Tabla comparativa de los intervalos de confianza y los parámetros reales
### obtenidos en el procedimiento Bootstrap
boot_comparison



