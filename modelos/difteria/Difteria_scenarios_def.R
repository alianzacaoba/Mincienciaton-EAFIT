#################################
# Libraries
#################################
library(deSolve)
library(dplyr)
library(ggplot2)
library(readxl)
options(scipen=999)
setwd("C:/Users/sulon/OneDrive - Universidad EAFIT/modelo-difteria/finales_1")

################################################################################
# Diphtheria SEIR+V model                                                      #
################################################################################
Difteria_30 = function(t, y, parms) {
  
  # Pull state variables from y vector
  S = y[1]
  I = y[2]
  R = y[3]
  V = y[4]
  W = y[5]
  
  # Pull parameter values from parms vector 
  ev1 = parms["ev1"] 
  gamma2 = parms["gamma2"] # tasa de recuperacion (sigma en paper)
  beta2 = parms["beta2"] # es k en el paper
  of = parms['of'] # 
  ii = parms['ii'] # 
  eta_v =parms['eta_v']
  eta_R =parms['eta_R']
  sd = parms['sd']
  vr = parms['vr']
  
  # social distancing flag
  dist <- (W>18330)*(W<18695)
  pdis <- (W>18330)*(W<(18330+sd*365))#marzo de 2020 al
  pvac <- (W>18330)*(W<(18330+sd*365))#marzo de 2020 al
  #dist <- 0
  
  # Fixed variables modeling (Based on population-birth_death_rates-vaccov-1960-2019 and meales_data_fit_models)
  # Birth rate
  brt <- ((W*-1.238913e-06 + 3.566146e-02))/365 # funcion de crecimiento Leandro, birth model
  # Death rate
  mu <- ((0.008251563  +  -4.65767e-07*W + -1.9377e-10 *W^2 + 7.939302e-14*W^3 + -1.160313e-17*W^4 + 8.236164e-22*W^5 + -2.837473e-26*W^6 + 3.794165e-31*W^7))/365
  # dose 1 vaccine coverage (base coverage input 1)
  cv1 <- 0.92*(1-exp(-0.07*(W)/365))*(W>11*365)*(1-vr*pvac)
  
  
  # modeling variability on the transmision rate
  betaD <- ( beta2*cos(2*pi*W/of) + beta2 )*(1-sd*pdis)
  
  # total population
  P <- S + I + R + V 
  
  # Define equations
  dS <- brt*(1-cv1*ev1)*P  - betaD*S*I/P - mu*S + eta_v*V + eta_R*R 
  dI <- betaD*S*I/P - (mu + gamma2)*I + ii
  dR <- gamma2 * I - mu*R - eta_R*R
  dV <- brt*cv1*ev1*P - eta_v*V - mu*V
  dW <- 1 
  dC <- (betaD*S*I)/P + ii 
  
  res = c(dS, dI, dR, dV, dW, dC)
  # Return list of gradients
  list(res)
}


################################################################################
#  SEIR Model parameters and initial conditions vectors for different sd and   #
#             vr scenarios.                                                    #
################################################################################
# This function receives two parameters:                                       #
# social_dist: amount of time social distancing was established in years       #
#               if you put 0.5 it'll simulate 6 months of social distancing    #
#                                                                              #
# vaccination_rate: amount decrease in vaccination rate during social          #
#                     distancing, if you put 0.1, it means vaccination rate    #
#                     decreased 10%                                            #
################################################################################

diphtheria_scenarios = function(social_dist, vaccination_rate) {
  
  # Parameters:
  prev1 <- 1 # Efectivity of dose 1
  prgamma2=1/14 # Tasa de recuperación 
  prbeta2 = 0.0594873493117294 # beta 
  # Initial state: y
  iniinfect = 226 # infectados iniciales (1997)
  pop = 21480065 # población de 1997
  iniI = iniinfect 
  iniR = 0
  iniS = iniI-iniR
  iniW = 0
  ii = 0.003
  
  parms = c(beta2 = prbeta2, gamma2 = prgamma2, 
            ev1=prev1, of = 365,ii=ii,eta_v=1/(365*30),
            eta_R=1/(365*30), vr = vaccination_rate, sd = social_dist) 
  start = c(S = pop, I = iniI, R = iniR, V=0, W = iniW, C=0)
  
  
  # time frame
  start.year <- 1970
  final.year <- 2029
  number.years <- final.year - start.year + 1
  end.time <- number.years*365
  time <- 1:end.time

  out <- ode(y=start, times=time, func=Difteria_30, 
             parms= parms)
  out <- as.data.frame(out) # model output as dataframe.
  
  
  casos_por_100k <-NULL # c?lcula base por a?os de casos por 100.000 habitantes y casos totales
  for (i in 1:number.years) { 
    fin = i*365
    #Casos
    tpop = out[fin,2]+out[fin,3]+out[fin,4]+out[fin,5]
    vac_a = out[fin,5]
    p_vac_a = vac_a/tpop
    cas = (out$C[i*365]-out$C[i*365-364])*100000/tpop
    cas0 = (out$C[i*365]-out$C[i*365-364])
    ## Totales
    rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0, tpop=tpop, vacunados=vac_a, prop_vac = p_vac_a)
    casos_por_100k<-rbind(casos_por_100k,rec)
  }
  
  return = c(out, casos_por_100k);
}

################################################################################
#  SEIR Model scenario simulation, with previous parameters and initial        #
#       conditions, varying social distancing and vaccination rate             #
################################################################################

#### Baseline: no social distancing, no decreased vaccination rate:

sd0_vr0 <- diphtheria_scenarios(0,0)

#### 6 months of social distancing - 
# 10% decrease in vaccination rate:
sd6_vr10 <- diphtheria_scenarios(0.5,0.1)

# 25% decrease in vaccination rate:
sd6_vr25 <- diphtheria_scenarios(0.5,0.25)

# 50% decrease in vaccination rate:
sd6_vr50 <- diphtheria_scenarios(0.5,0.5)


#### 12 months of social distancing - 
# 10% decrease in vaccination rate:
sd12_vr10 <- diphtheria_scenarios(1,0.1)

# 25% decrease in vaccination rate:
sd12_vr25 <- diphtheria_scenarios(1,0.25)

# 50% decrease in vaccination rate:
sd12_vr50 <- diphtheria_scenarios(1,0.5)


################################################################################
#  Scenarios visualization                                                     #
################################################################################

####  Anual cases for all simulated period (2019-2029):
par(mfrow=c(1,2))

plot(2019:2029,sd0_vr0$casos_t[50:60],type="l",lty=1,col="black",lwd=2,xlab="",
     ylab="Annual cases",main="Diphtheria anual cases (2019-2029)",ylim=c(0,max(sd0_vr0$casos_t[50:60])),las=1,xaxt="none")
axis(1, seq(2019,2029,1))
lines(2019:2029,sd6_vr10$casos_t[50:60],type="l",col="blue",lwd=2)
lines(2019:2029,sd6_vr25$casos_t[50:60],type="l",col="green",lwd=2)
lines(2019:2029,sd6_vr50$casos_t[50:60],type="l",col="red",lwd=2)
lines(2019:2029,sd12_vr10$casos_t[50:60],type="l",col="blue",lwd=2,lty=2)
lines(2019:2029,sd12_vr25$casos_t[50:60],type="l",col="green",lwd=2,lty=2)
lines(2019:2029,sd12_vr50$casos_t[50:60],type="l",col="red",lwd=2,lty=2)
legend("bottomright", inset=c(-0.85,0), legend=c("Baseline","Six months (coverage -10%)","Six months (coverage -25%)","Six months (coverage -50%)",
                                                 "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.9,  bty = 'n')

#### Cases per 100k inhabitants from 2024 - 2029 

plot(2024:2029,sd0_vr0$casos_t[55:60],type="l",lty=1,col="black",lwd=2,xlab="",
     ylab="Anual cases",main="Diphtheria anual cases (2024-2029)",
     ylim=c(min(sd0_vr0$casos_t[55:60]),max(sd12_vr50$casos_t[55:60])),las=1)
lines(2024:2029,sd6_vr10$casos_t[55:60],type="l",col="blue",lwd=2)
lines(2024:2029,sd6_vr25$casos_t[55:60],type="l",col="green",lwd=2)
lines(2024:2029,sd6_vr50$casos_t[55:60],type="l",col="red",lwd=2)
lines(2024:2029,sd12_vr10$casos_t[55:60],type="l",col="blue",lwd=2,lty=2)
lines(2024:2029,sd12_vr25$casos_t[55:60],type="l",col="green",lwd=2,lty=2)
lines(2024:2029,sd12_vr50$casos_t[55:60],type="l",col="red",lwd=2,lty=2)
# legend("bottomright", inset=c(-0.85,0), legend=c("Baseline","Six months (coverage -10%)","Six months (coverage -25%)","Six months (coverage -50%)",
#                                               "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
#        col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.9, bty = 'n')


################################################################################
#                               Con protecci?n                                 #
################################################################################

anualizar_prop_vaunas <- function(prop_vacunas) {
  PROP_VAC_ANUAL<-data.frame(ano = numeric(0),prop_vac_anual = numeric(0)) # cálcula base por años de casos por 100.000 habitantes y casos totales
  for (i in 1:60) { 
    fin = i*365
    #vac_anual = data.frame(ano=i+1969, prop_vac_anual=sum(prop_vacunas[(((i-1)*365)+1):(i*365)]))
    vac_anual = data.frame(ano=i+1969, prop_vac_anual=sum(prop_vacunas[i*365]))
    PROP_VAC_ANUAL = rbind(PROP_VAC_ANUAL, vac_anual)
    
  }
  return(PROP_VAC_ANUAL)
}

# Total personas en el modelo

##### Proporci?n de protegidos sin COVID-19
P01 <- sd0_vr0$V/(sd0_vr0$S + sd0_vr0$I + sd0_vr0$R + sd0_vr0$V) # MMR1

##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 10%
P6_101<- sd6_vr10$V/(sd6_vr10$S + sd6_vr10$I + sd6_vr10$R + sd6_vr10$V)# MMR1
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 25%
P6_251<- sd6_vr25$V/(sd6_vr25$S + sd6_vr25$I + sd6_vr25$R + sd6_vr25$V) # MMR1
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 50%
P6_501<- sd6_vr50$V/(sd6_vr50$S + sd6_vr50$I + sd6_vr50$R + sd6_vr50$V) # MMR1

##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 10%
P1_101<- sd12_vr10$V/(sd12_vr10$S + sd12_vr10$I + sd12_vr10$R + sd12_vr10$V) # MMR1
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 25%
P1_251<- sd12_vr25$V/(sd12_vr25$S + sd12_vr25$I + sd12_vr25$R + sd12_vr25$V) # MMR1
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 50%
P1_501<- sd12_vr50$V/(sd12_vr50$S + sd12_vr50$I + sd12_vr50$R + sd12_vr50$V) # MMR1


##### Proporci?n de protegidos sin COVID-19
P01_a <- anualizar_prop_vaunas(P01) # DTP

##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 10%
P6_101_a <- anualizar_prop_vaunas(P6_101) # DTP
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 25%
P6_251_a <- anualizar_prop_vaunas(P6_251) # DTP
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 50%
P6_501_a <- anualizar_prop_vaunas(P6_501) # DTP

##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 10%
P1_101_a <- anualizar_prop_vaunas(P1_101) # DTP
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 25%
P1_251_a <- anualizar_prop_vaunas(P1_251) # DTP
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 50%
P1_501_a <- anualizar_prop_vaunas(P1_501) # DTP


##################################### Figura proporci?n de protegidos

#### Tiempo en d?as entre el 2019 y 2029
Te=18250:21900
# par(mfrow=c(1,1))
# plot(P01_a$ano[50:60],P01_a$prop_vac_anual[50:60],type="l",lty=1,col="black",lwd=2,
#      ylim=c(min(P1_501_a$prop_vac_anual[50:60]),max(P01_a$prop_vac_anual[50:60])),xlab ="",ylab="Proportion",
#      main="Protected by DTP 3 (Diphtheria)") 
# lines(P1_101_a$ano[50:60],P1_101_a$prop_vac_anual[50:60],type="l",col="blue",lwd=2)
# lines(P1_251_a$ano[50:60],P1_251_a$prop_vac_anual[50:60],type="l",col="green",lwd=2)
# lines(P1_501_a$ano[50:60],P1_501_a$prop_vac_anual[50:60],type="l",col="red",lwd=2)
# lines(P6_101_a$ano[50:60],P6_101_a$prop_vac_anual[50:60],type="l",col="blue",lwd=2,lty=2)
# lines(P6_251_a$ano[50:60],P6_251_a$prop_vac_anual[50:60],type="l",col="green",lwd=2,lty=2)
# lines(P6_501_a$ano[50:60],P6_501_a$prop_vac_anual[50:60],type="l",col="red",lwd=2,lty=2)
par(mfrow=c(1,1))
plot(Te/365+1969,P01[18250:21900],type="l",lty=1,col="black",lwd=2,
     ylim=c(min(P1_501[18250:21900]),max(P01[18250:21900])),xlab ="",ylab="Proportion",
     main="Protected by DPT3 (Diphtheria)")
lines(Te/365+1969,P6_101[18250:21900],type="l",col="blue",lwd=2)
lines(Te/365+1969,P6_251[18250:21900],type="l",col="green",lwd=2)
lines(Te/365+1969,P6_501[18250:21900],type="l",col="red",lwd=2)
lines(Te/365+1969,P1_101[18250:21900],type="l",col="blue",lwd=2,lty=2)
lines(Te/365+1969,P1_251[18250:21900],type="l",col="green",lwd=2,lty=2)
lines(Te/365+1969,P1_501[18250:21900],type="l",col="red",lwd=2,lty=2)

legend("bottomleft", legend=c("Baseline","Six months (coverage -10%)","Six months (coverage -25%)","Six months (coverage -50%)",
                            "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.8, bty = 'n')


####################################################################################################
#                                                                                                  #
#                               Datos de salida del modelo                                         #
#                                                                                                  #
####################################################################################################
# Pronóstico casos
SC<-data.frame("id-variable"=rep("UE-DIF-SC",11),"Unidades temporales"=sd0_vr0$ano[50:60],Valor_proyectado=sd0_vr0$casos_t[50:60])
CV610<-data.frame("id-variable"=rep("UE-DIF-CV610",11),"Unidades temporales"=sd6_vr10$ano[50:60],Valor_proyectado=sd6_vr10$casos_t[50:60])
CV625<-data.frame("id-variable"=rep("UE-DIF-CV625",11),"Unidades temporales"=sd6_vr25$ano[50:60],Valor_proyectado=sd6_vr25$casos_t[50:60])
CV650<-data.frame("id-variable"=rep("UE-DIF-CV650",11),"Unidades temporales"=sd6_vr50$ano[50:60],Valor_proyectado=sd6_vr50$casos_t[50:60])
CV110<-data.frame("id-variable"=rep("UE-DIF-CV110",11),"Unidades temporales"=sd12_vr10$ano[50:60],Valor_proyectado=sd12_vr10$casos_t[50:60])
CV125<-data.frame("id-variable"=rep("UE-DIF-CV125",11),"Unidades temporales"=sd12_vr10$ano[50:60],Valor_proyectado=sd12_vr10$casos_t[50:60])
CV150<-data.frame("id-variable"=rep("UE-DIF-CV150",11),"Unidades temporales"=sd12_vr10$ano[50:60],Valor_proyectado=sd12_vr10$casos_t[50:60])
BASE_PRO<-rbind(SC,CV610,CV625,CV650,CV110,CV125,CV150)
library("xlsx")
#write.xlsx(BASE_PRO, "RUB_Pronostico_scenarios.xlsx")

# Pronóstico proporción de protegidos para la rubeola (vacuna MMR1) 
V1SC<-data.frame("id-variable"=rep("UE-DIF-PV1SC",11),"Unidades temporales"=2019:2029,Valor_proyectado=P01_a[50:60,])
V1610<-data.frame("id-variable"=rep("UE-DIF-PV1610",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_101_a[50:60,])
V1625<-data.frame("id-variable"=rep("UE-DIF-PV1625",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_251_a[50:60,])
V1650<-data.frame("id-variable"=rep("UE-DIF-PV1650",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_501_a[50:60,])
V1110<-data.frame("id-variable"=rep("UE-DIF-PV1110",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_101_a[50:60,])
V1125<-data.frame("id-variable"=rep("UE-DIF-PV1125",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_251_a[50:60,])
V1150<-data.frame("id-variable"=rep("UE-DIF-PV1150",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_501_a[50:60,])
