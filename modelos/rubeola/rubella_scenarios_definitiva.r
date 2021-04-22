################################################################################
#                     Rubella SEIR + 2 dose vaccination                        #
# this script is a SEIRV2 model and a function to change social distancing and #
# vaccination rate.                                                            #
# It was created to understand the impact of COVID-19 impact on Rubella cases  #
# and if the decreased vaccination coverage has an impact.                     #  
################################################################################

################################################################################
#  Libraries                                                                   #
################################################################################

library(xlsx)
library(deSolve)
library(dfoptim)
options(scipen=999)

################################################################################
#  SEIR Model function                                                         #
################################################################################

#### Model function:
# t is the time frame
# y is the initial conditions vector 
# parms are the model parameters 

#### time frame
# start.year <- 1970 # initial year 
# final.year <- 2021 # last year

#### Initial conditions:
# iniS: initial Suceptibles.
# iniE: initial Exposen.
# iniI: initial Infected.
# iniR: initial Recovered.
# iniV: initial dose 1 vaccinated.
# iniV2: initial dose 2 vaccinated.
# iniW: initial time.
# iniC: initial cases. 

# initial_cond_vector <- c(S = iniS, E = iniE, I = iniI, R = iniR, 
#                          V=0, V2=0, W = iniW, C=iniC)

#### Needed parameters for the model:
# beta2: rate at which Suceptibles become Exposed.
# gamma2: rate at which infected become Recovered.
# ev1: effectivity of vaccine, dose 1.
# ev2: effectivity of vaccine, dose 2.
# ee: rate at which newborns of Exposed born into Exposed.
# ie: rate at which newborns of Exposed born into Infected.
# e: rate at which exposed become infected.
# of: period at which new infected are entered into the model. 
# ii: imported cases per year. 
# sd:  
# vr:  

# parameters_vector <- c(beta2 = prbeta2, gamma2 = prgamma2, ev1=prev1, 
#                        ev2 = prev2, ee = pree, ie = prie,
#                        e = pre, of=prof, ii=prii)

#### For running the model:
# identifier <- ode(y = initial_cond_vector, 
#                   times = time_frame_vector, 
#                   func = model_definition_function, 
#                   parms = parameter_vector)


# Model function:
rubella_model <- function(t, y, parms) {
  # Pull initial conditions y vector
  S = y[1]
  E = y[2]
  I = y[3]
  R = y[4]
  V = y[5]
  V2 = y[6]
  W = y[7]
  
  # Pull parameter values from parms vector 
  ev1 = parms["ev1"] 
  ev2 = parms["ev2"]
  gamma2 = parms["gamma2"] # tasa de recuperacion (sigma en paper)
  N = parms["N"] # total poblacion
  beta2 = parms["beta2"] # es k en el paper
  ee = parms["ee"] # en paper es p (hijos de los de la clase E que entran a E)
  ie = parms["ie"]  # en paper es q (hijos de los de la clase I que entran a E)
  e = parms["e"] # tasa a la cual los expuestos se vuelven infectados
  of = parms['of'] # 
  ii = parms['ii'] # 
  sd = parms['sd']
  vr = parms['vr']
  
  #### social distancing flag
  # dist <- (W>18330)*(W<18695)
  # dist <- 0
  pdis = (W>18330)*(W<(18330+sd*365))#marzo de 2020 al 
  pvac = (W>18330)*(W<(18330+sd*365))#marzo de 2020 al
  
  
  
  # Fixed variables modeling (Based on 
  # population-birth_death_rates-vaccov-1960-2019 and meales_data_fit_models)
  #### Birth rate
  brt <- ((W*-1.238913e-06 + 3.566146e-02))/365 
  
  #### Death rate
  mu <- ((0.008251563  +  -4.65767e-07*W + -1.9377e-10 *W^2 + 
            7.939302e-14*W^3 + -1.160313e-17*W^4 + 8.236164e-22*W^5 + 
            (-2.837473e-26)*W^6 + 3.794165e-31*W^7))/365
  
  #### dose 1 vaccine coverage (base coverage input 1)
  cv1 <- (W < 3651)*0 + (W >= 3651)*(0.95*(1-exp(-0.1408442*(W-3284)/365)))*(1-vr*pvac) #+ 
    #(0.9-vr)*(W>=18695) + 0.9*(W>18330)*(W<18695)
  #cv1=0
  
  #### dose 2 vaccine coverage
  cv2 <- (W < 9856)*0 + (W >= 9856)*(0.85-exp(-0.4*(W-9489)/365))*(1-vr*pvac) #+
    #(0.9-vr)*(W>=18695) + 0.9*(W>18330)*(W<18695)
  #cv2=0
  
  #### modeling variability on the transmision rate
  betaD <- ( beta2*cos(2*pi*W/of) + beta2 )*(1-0.4*pdis)
  
  #### total population
  P <- S + E + I + R + V + V2
  
  # Define equations
  dS <- brt*(1-cv1*ev1)*P + V*(1/(365*4))*cv2*(1-ev2) - betaD*S*I/P - ee*brt*E - ie*brt*E - mu*S #br*N - beta*S*I/N - br*N*v1c - mr*S
  dE <- betaD*S*I/P + ee*brt*E + ie*brt*E - e*E - mu*E
  dI <- e*E - (mu + gamma2)*I + ii*(W>=18616) + 0*(W>=18616)*(W<=18617) 
  dR <- gamma2 * I - mu*R
  dV <- brt*cv1*ev1*P - V*(1/(365*4))*cv2 - mu*V
  dV2 <- V*(1/(365*4))*cv2*ev2 - mu*V2
  dW <- 1 
  dC <- e*E + ii*(W>=18616) + 0*(W>=18616)*(W<=18617) # new cases
  
  res = c(dS, dE, dI, dR, dV, dV2, dW, dC)
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



rubella_scenarios = function(social_dist, vaccination_rate) {
  
  #### time frame
  start.year <- 1970
  final.year <- 2029
  number.years <- final.year - start.year + 1
  end.time <- number.years*365
  times <- 1:end.time
  
  
  #### Needed parameters for the model:
  # beta2: rate at which Suceptibles become Exposed.
  # gamma2: rate at which infected become Recovered.
  # ev1: effectivity of vaccine, dose 1.
  # ev2: effectivity of vaccine, dose 2.
  # ee: rate at which newborns of Exposed born into Exposed.
  # ie: rate at which newborns of Exposed born into Infected.
  # e: rate at which exposed become infected.
  # of: period at which new infected are entered into the model. 
  # ii: imported cases per year. 
  # sd: social distancing in years.
  # vr: vaccination rate decrease
  
  # parms = c(beta2 = prbeta2, gamma2 = prgamma2, ev1=prev1, 
  #           ev2 = prev2, ee = pree, ie = prie,
  #           e = pre, of=prof, ii=prii, sd = social_dist, 
  #           vr = vaccination_rate)
  # beta_optimo = 0.01681578 
  # gamma_optimo = 0.01340360
  
  prbeta2  <- 0.01681578
  prgamma2 <- 0.01340360
  prev1    <- 0.88
  prev2    <- 0.97
  pree     <- 0.65 # de acuerdo al paper
  prie     <- 0.65 # de acuerdo al paper
  pre      <-  1/90 # de acuerdo al paper
  prof     <- 365
  prii     <- 0/365
  
  
  parms = c(beta2 = prbeta2, gamma2 = prgamma2, 
            ev1=prev1, ev2 = prev2, ee = pree, 
            ie = prie, e = pre, of=prof, ii=prii,
            sd = social_dist, vr = vaccination_rate)
  
  #### Initial conditions:
  # iniS: initial Suceptibles.
  # iniE: initial Exposen.
  # iniI: initial Infected.
  # iniR: initial Recovered.
  # iniV: initial dose 1 vaccinated.
  # iniV2: initial dose 2 vaccinated.
  # iniW: initial time.
  # iniC: initial cases. 
  
  # start = c(S = iniS, E = iniE, I = iniI, R = iniR, 
  #           V=0, V2=0, W = iniW, C=iniC)
  
  iniinfect = 5 # initial number of exposed.
  pop = 21480065 # population in Colombia in 1970
  iniI = iniinfect 
  iniR = 0
  iniE = 0
  iniS = pop-iniI-iniR
  iniW = 1
  iniC = 0
  iniV = 0
  iniV2 = 0
  
  start = c(S = iniS, E = iniE, I = iniI, R = iniR, 
            V=iniV, V2=iniV2, W = iniW, C=iniC)
  
  #### For running the model:
  # identifier <- ode(y = initial_cond_vector, 
  #                   times = time_frame_vector, 
  #                   func = model_definition_function, 
  #                   parms = parameter_vector)
  
  out <- ode(y=start, times=times, func=rubella_model, parms= parms)
  out <- as.data.frame(out) # model output as dataframe.
  
  casos_por_100k <-NULL # c?lcula base por a?os de casos por 100.000 habitantes y casos totales
  for (i in 1:60) { 
    fin = i*365
    #Casos
    tpop = out[fin,2]+out[fin,3]+out[fin,4]+out[fin,5]+out[fin,6]+out[fin,7]
    cas = (out$C[i*365]-out$C[i*365-364])*100000/tpop
    cas0 = (out$C[i*365]-out$C[i*365-364])
    ## Totales
    rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0)
    casos_por_100k<-rbind(casos_por_100k,rec)
  }
  
  return = c(out, casos_por_100k);
}

################################################################################
#  SEIR Model scenario simulation, with previous parameters and initial        #
#       conditions, varying social distancing and vaccination rate             #
################################################################################

#### Baseline: no social distancing, no decreased vaccination rate:
sd0_vr0 <- rubella_scenarios(0,0)

#### 6 months of social distancing - 
# 10% decrease in vaccination rate:
sd6_vr10 <- rubella_scenarios(0.5,0.1)

# 25% decrease in vaccination rate:
sd6_vr25 <- rubella_scenarios(0.5,0.25)

# 50% decrease in vaccination rate:
sd6_vr50 <- rubella_scenarios(0.5,0.5)


#### 12 months of social distancing - 
# 10% decrease in vaccination rate:
sd12_vr10 <- rubella_scenarios(1,0.1)

# 25% decrease in vaccination rate:
sd12_vr25 <- rubella_scenarios(1,0.25)

# 50% decrease in vaccination rate:
sd12_vr50 <- rubella_scenarios(1,0.5)


################################################################################
#  Scenarios visualization                                                     #
################################################################################

####  Anual cases for all simulated period (2019-2029):
par(mfrow=c(1,2))

plot(2019:2029,sd0_vr0$casos_t[50:60],type="l",lty=1,col="black",lwd=2,xlab="",
     ylab="Annual cases",main="Rubella anual cases (2019-2029)",ylim=c(0,max(sd0_vr0$casos_t[50:60])),las=1,xaxt="none")
axis(1, seq(2019,2029,1))
lines(2019:2029,sd6_vr10$casos_t[50:60],type="l",col="blue",lwd=2)
lines(2019:2029,sd6_vr25$casos_t[50:60],type="l",col="green",lwd=2)
lines(2019:2029,sd6_vr50$casos_t[50:60],type="l",col="red",lwd=2)
lines(2019:2029,sd12_vr10$casos_t[50:60],type="l",col="blue",lwd=2,lty=2)
lines(2019:2029,sd12_vr25$casos_t[50:60],type="l",col="green",lwd=2,lty=2)
lines(2019:2029,sd12_vr50$casos_t[50:60],type="l",col="red",lwd=2,lty=2)
#legend("topright", inset=c(-0.2,0), legend=c("Baseline","Six months (coverage -10%)","Six months (coverage -25%)","Six months (coverage -50%)",
#                             "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
#       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.9)

#### Cases per 100k inhabitants from 2024 - 2029 

plot(2024:2029,sd0_vr0$casos_t[55:60],type="l",lty=1,col="black",lwd=2,xlab="",
     ylab="Anual cases",main="Rubella anual cases (2024-2029)",ylim=c(min(sd0_vr0$casos_t[55:60]),
                                                                      max(sd0_vr0$casos_t[55:60])),las=1)
lines(2024:2029,sd6_vr10$casos_t[55:60],type="l",col="blue",lwd=2)
lines(2024:2029,sd6_vr25$casos_t[55:60],type="l",col="green",lwd=2)
lines(2024:2029,sd6_vr50$casos_t[55:60],type="l",col="red",lwd=2)
lines(2024:2029,sd12_vr10$casos_t[55:60],type="l",col="blue",lwd=2,lty=2)
lines(2024:2029,sd12_vr25$casos_t[55:60],type="l",col="green",lwd=2,lty=2)
lines(2024:2029,sd12_vr50$casos_t[55:60],type="l",col="red",lwd=2,lty=2)
legend("topright", inset=c(-0.85,0), legend=c("Baseline","Six months (coverage -10%)","Six months (coverage -25%)","Six months (coverage -50%)",
                                              "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.9, bty = 'n')

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
Total <- (sd0_vr0$S+sd0_vr0$E+sd0_vr0$I+sd0_vr0$R+
            sd0_vr0$V+sd0_vr0$V2)

##### Proporci?n de protegidos sin COVID-19
P01 <- sd0_vr0$V/Total # MMR1
P02 <- sd0_vr0$V2/Total # MMR2

##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 10%
P6_101<- sd6_vr10$V/Total # MMR1
P6_102<- sd6_vr10$V2/Total # MMR2
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 25%
P6_251<- sd6_vr25$V/Total # MMR1
P6_252<- sd6_vr25$V2/Total # MMR2
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 50%
P6_501<- sd6_vr50$V/Total # MMR1
P6_502<- sd6_vr50$V2/Total # MMR2

##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 10%
P1_101<- sd12_vr10$V/Total # MMR1
P1_102<- sd12_vr10$V2/Total # MMR2
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 25%
P1_251<- sd12_vr25$V/Total # MMR1
P1_252<- sd12_vr25$V2/Total # MMR2
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 50%
P1_501<- sd12_vr50$V/Total # MMR1
P1_502<- sd12_vr50$V2/Total # MMR2



##### Proporci?n de protegidos sin COVID-19
P01_a <- anualizar_prop_vaunas(sd0_vr0$V/Total) # MMR1
P02_a <- anualizar_prop_vaunas(sd0_vr0$V2/Total) # MMR2

##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 10%
P6_101_a <- anualizar_prop_vaunas(sd6_vr10$V/Total) # MMR1
P6_102_a <-anualizar_prop_vaunas(sd6_vr10$V2/Total) # MMR2
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 25%
P6_251_a <- anualizar_prop_vaunas(sd6_vr25$V/Total) # MMR1
P6_252_a <- anualizar_prop_vaunas(sd6_vr25$V2/Total) # MMR2
##### Proporci?n de protegidos con 6 meses de distanciamiento con caida en cobertura del 50%
P6_501_a <- anualizar_prop_vaunas(sd6_vr50$V/Total) # MMR1
P6_502_a <- anualizar_prop_vaunas(sd6_vr50$V2/Total) # MMR2

##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 10%
P1_101_a <- anualizar_prop_vaunas(sd12_vr10$V/Total) # MMR1
P1_102_a <- anualizar_prop_vaunas(sd12_vr10$V2/Total) # MMR2
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 25%
P1_251_a <- anualizar_prop_vaunas(sd12_vr25$V/Total) # MMR1
P1_252_a <- anualizar_prop_vaunas(sd12_vr25$V2/Total) # MMR2
##### Proporci?n de protegidos con un a?o de distanciamiento con caida en cobertura del 50%
P1_501_a <- anualizar_prop_vaunas(sd12_vr50$V/Total) # MMR1
P1_502_a <- anualizar_prop_vaunas(sd12_vr50$V2/Total) # MMR2

#### Tiempo en d?as entre el 2019 y 2029
Te=18250:21900


##################################### Figura proporci?n de protegidos

par(mfrow=c(1,2))

plot(Te/365+1969,P01[18250:21900],type="l",lty=1,col="black",lwd=2,
     ylim=c(min(P1_501[18250:21900]),max(P01[18250:21900])),xlab ="",ylab="Proportion",
     main="Protected by MMR 1 (rubella)")
lines(Te/365+1969,P6_101[18250:21900],type="l",col="blue",lwd=2)
lines(Te/365+1969,P6_251[18250:21900],type="l",col="green",lwd=2)
lines(Te/365+1969,P6_501[18250:21900],type="l",col="red",lwd=2)
lines(Te/365+1969,P1_101[18250:21900],type="l",col="blue",lwd=2,lty=2)
lines(Te/365+1969,P1_251[18250:21900],type="l",col="green",lwd=2,lty=2)
lines(Te/365+1969,P1_501[18250:21900],type="l",col="red",lwd=2,lty=2)

#legend("topright", legend=c("Baseline","Six months (coverage -10%)","Six months (coverage -25%)","Six months (coverage -50%)",
#                            "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
#       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.9)


plot(Te/365+1969,P02[18250:21900],type="l",lty=1,col="black",lwd=2,
     ylim=c(min(P1_502[18250:21900]),max(P02[18250:21900])),xlab ="",
     ylab="Proportion",main="Protected by MMR 2 (rubella)")
lines(Te/365+1969,P6_102[18250:21900],type="l",col="blue",lwd=2)
lines(Te/365+1969,P6_252[18250:21900],type="l",col="green",lwd=2)
lines(Te/365+1969,P6_502[18250:21900],type="l",col="red",lwd=2)
lines(Te/365+1969,P1_102[18250:21900],type="l",col="blue",lwd=2,lty=2)
lines(Te/365+1969,P1_252[18250:21900],type="l",col="green",lwd=2,lty=2)
lines(Te/365+1969,P1_502[18250:21900],type="l",col="red",lwd=2,lty=2)
legend("topleft", legend=c("Baseline","Six months (coverage -10%)","Six months (coverage -25%)","Six months (coverage -50%)",
                            "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.8, bty = 'n')


####################################################################################################
#                                                                                                  #
#                               Datos de salida del modelo                                         #
#                                                                                                  #
####################################################################################################
# Pronóstico casos
SC<-data.frame("id-variable"=rep("UE-RUB-SC",11),"Unidades temporales"=sd0_vr0$ano[50:60],Valor_proyectado=sd0_vr0$casos_t[50:60])
CV610<-data.frame("id-variable"=rep("UE-RUB-CV610",11),"Unidades temporales"=sd6_vr10$ano[50:60],Valor_proyectado=sd6_vr10$casos_t[50:60])
CV625<-data.frame("id-variable"=rep("UE-RUB-CV625",11),"Unidades temporales"=sd6_vr25$ano[50:60],Valor_proyectado=sd6_vr25$casos_t[50:60])
CV650<-data.frame("id-variable"=rep("UE-RUB-CV650",11),"Unidades temporales"=sd6_vr50$ano[50:60],Valor_proyectado=sd6_vr50$casos_t[50:60])
CV110<-data.frame("id-variable"=rep("UE-RUB-CV110",11),"Unidades temporales"=sd12_vr10$ano[50:60],Valor_proyectado=sd12_vr10$casos_t[50:60])
CV125<-data.frame("id-variable"=rep("UE-RUB-CV125",11),"Unidades temporales"=sd12_vr10$ano[50:60],Valor_proyectado=sd12_vr10$casos_t[50:60])
CV150<-data.frame("id-variable"=rep("UE-RUB-CV150",11),"Unidades temporales"=sd12_vr10$ano[50:60],Valor_proyectado=sd12_vr10$casos_t[50:60])
BASE_PRO<-rbind(SC,CV610,CV625,CV650,CV110,CV125,CV150)
library("xlsx")
write.xlsx(BASE_PRO, "RUB_Pronostico_scenarios.xlsx")

# Pronóstico proporción de protegidos para la rubeola (vacuna MMR1) 
V1SC<-data.frame("id-variable"=rep("UE-RUB-PV1SC",11),"Unidades temporales"=2019:2029,Valor_proyectado=P01_a[50:60,])
V1610<-data.frame("id-variable"=rep("UE-RUB-PV1610",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_101_a[50:60,])
V1625<-data.frame("id-variable"=rep("UE-RUB-PV1625",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_251_a[50:60,])
V1650<-data.frame("id-variable"=rep("UE-RUB-PV1650",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_501_a[50:60,])
V1110<-data.frame("id-variable"=rep("UE-RUB-PV1110",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_101_a[50:60,])
V1125<-data.frame("id-variable"=rep("UE-RUB-PV1125",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_251_a[50:60,])
V1150<-data.frame("id-variable"=rep("UE-RUB-PV1150",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_501_a[50:60,])



# Pronóstico proporción de protegidos para la rubeola (vacuna MMR2) 
V2SC<-data.frame("id-variable"=rep("UE-RUB-PV2SC",11),"Unidades temporales"=2019:2029,Valor_proyectado=P02_a[50:60,])
V2610<-data.frame("id-variable"=rep("UE-RUB-PV2610",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_102_a[50:60,])
V2625<-data.frame("id-variable"=rep("UE-RUB-PV2625",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_252_a[50:60,])
V2650<-data.frame("id-variable"=rep("UE-RUB-PV2650",11),"Unidades temporales"=2019:2029,Valor_proyectado=P6_502_a[50:60,])
V2110<-data.frame("id-variable"=rep("UE-RUB-PV2110",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_102_a[50:60,])
V2125<-data.frame("id-variable"=rep("UE-RUB-PV2125",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_252_a[50:60,])
V2150<-data.frame("id-variable"=rep("UE-RUB-PV2150",11),"Unidades temporales"=2019:2029,Valor_proyectado=P1_502_a[50:60,])
# Base completa de protegidos
BASE_VAC<-rbind(V1SC,V1610,V1625,V1650,V1110,V1125,V1150,V2SC,V2610,V2625,V2650,V2110,V2125,V2150)
write.xlsx(BASE_VAC, "RUB_Pronostico_scenarios_Vacc.xlsx")
