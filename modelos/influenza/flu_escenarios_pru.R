### Escenarios de simulacion
### Funcion
library(deSolve)
# Modelo
flu3 <- function(time, y, parms) {
  with(as.list(c(y, parms)),{
    
    ###############
    pdis = (W>585)*(W<(585+tdis*52))#marzo de 2020 al 
    pvac = (W>585)*(W<(585+tdis*52))#marzo de 2020 al 
    
    W2=W
    mu2 = ((4.347058e-01 -5.035106e-03*W2 + 2.982051e-05*W2^2 + 7.823542e-07*W2^3 -8.044537e-09*W2^4 + 2.997302e-11*W2^5 -5.035786e-14*W2^6 + 3.212231e-17*W2^7)*(W2<(9*52)) + 0.7058571*(W2>(9*52-1)))*(1-rvac*pvac*va) # Funcion de cobertura del 2009 en adelante
    mu3 = (0.6166667*(W2>(52*4)))*(1-rvac*pdis*va) # Cobertura de vacunacion de grupo de riesgo 2 (gestantes) del 2009 en adelante
    mu4 = 0.5*(1-rvac*pvac*va)#((2.225411e-01 +  3.111241e-02*W2 -6.578987e-04*W2^2 + 4.361894e-06*W2^3 -9.084438e-09*W2^4 -3.937422e-12*W2^5 + 2.165713e-14*W2^6)*(W2<(7*52)) + 0.71*(W2>(7*52-1)))*(1-rvac*pdis*va)# cobertura de vacunacion grupo de riesgo 3 (adultos mayores)
    
    B1 = (B*cos(2*pi/52*(W-of))+B)*(1-0.4*pdis*dist)# tasa de trasmision de susceptibles totales
    N <- S1+V1P+V1G+E1+I1+R1+S2+V2+E2+I2+R2 # Número total de personas
    
    ############# Entrada de casos externos migracion despues del 2015 (0 migrantes por 100 mil) que es positiva hasta llegar a 4.16 migrantes por 100 mil en 2020
    la2<-im2*((W*0.016-4.992)/1000)/52*N*(W>312)*(W<624)*(1-0.4*pdis*dist)+im2*((4.16)/1000)/52*N*(W>624)*(W<728)*(1-0.4*pdis*dist)
    
    ##### Crecimiento poblacional
    
    n <- ((W1*-1.238913e-06 + 3.566146e-02))/52 #funcion de crecimiento poblacional
    mu <-((0.008251563  +  -4.65767e-07*W1 + -1.9377e-10 *W1^2 + 7.939302e-14*W1^3 + -1.160313e-17*W1^4 + 8.236164e-22*W1^5 + -2.837473e-26*W1^6 + 3.794165e-31*W1^7))/52 #funcion de tasa de muerte
    
    # Estados
    
    dS1 <- N*n*(1-mu2*ef1) +rho*R1 +tau*V1P +tau*V1G -B1*S1*(I1+I2)/N - N*n*mu3*ef2 - phi*S1 - mu*S1 - im*(1-0.4*pdis*dist) # susceptibles antes de los 60
    dV1P <- N*n*mu2*ef1 -tau*V1P -mu*V1P # vacunados grupo 1 (ninos menores a un ano)
    dV1G <- N*n*mu3*ef2 -tau*V1G -mu*V1G # vacunados grupo 2 (madres gestantes)
    dE1 <-  B1*S1*(I1+I2)/N - al*E1 -mu*E1 +im*(1-0.4*pdis*dist) -phi*E1# Expuestos menores de 60 anos
    dI1 <-  al*E1 -ga*I1 -mu*I1 -mud*I1 +la2*(1-0.066)*(1-pdis*dist) -phi*I1# Infectados menores de 60 anos
    dR1 <-  ga*I1 -rho*R1 -mu*R1 -phi*R1# Recuperados menores de 60 anos
    dS2 <- phi*S1 +rho*R2 +tau*V2 -B1*S2*(I1+I2)/N - mu4*ef3*S2 - mu*S2 -im*(1-0.4*pdis*dist)# susceptibles mayores de 60
    dV2 <- mu4*ef3*S2 -tau*V2 -mu*V2 # vacunados grupo 3 (mayores de 60 anos)
    dE2 <- B1*S2*(I1+I2)/N -al*E2 -mu*E2 +im*(1-0.4*pdis*dist) +phi*E1# Expuestos menores de 60 anos
    dI2 <- al*E2 -ga*I2 -mu*I2 -mud*I2 + la2*0.066*(1-pdis*dist) +phi*I1 # Infectados menores de 60 anos
    dR2 <- ga*I2 -rho*R2 -mu*R2 -phi*R1# Recuperados menores de 60 anos
    
    ############## Adicionales para calculos
    
    dW = 1 # inicio en 2009  # tiempo de simulacion # temporizador para funciones de vacunacion y de tasa de transmision
    dW1 = 1 # temporizador para funciones de tasa de nacimiento y de mortalidad
    dC = al*(E1+E2) +la2*(1-pdis*dist) # Nuevos casos # transmission (al*E) e importados (im)
    dCV1 = N*n*mu2*ef1
    dCV2 = N*n*mu3*ef1
    dCV3 = mu4*ef1*S2
    
    list(c(dS1,dV1P,dV1G,dE1,dI1,dR1,dS2,dV2,dE2,dI2,dR2,dW,dW1,dC,dCV1,dCV2,dCV3))})}

######################################
## Funcion para simular escenarios   #
######################################
# Pide tiempo de distanciamiento (tdis)
# y porcentaje de caida en la cobertura
# de vacunacion durante la pandemia (rvac)

simuflu = function(tdis,rvac) {
  
  parms <- c(
    B=0.732169990,#0.15*7,#
    ef1 = 0.74, # efectividad vacuna pop1
    ef2 = 0.54, # efectividad vacuna pop2
    ef3 = 0.4, # efectividad vacuna pop3
    tau = 1/52, # tasa de perdida de inmunidad vacuna 
    al = 1/(1/7),# tasa de paso de expuestos a infectados
    ga = 1/(5/7), # tasa de paso de infectados a recuperados
    im = 1.905911899,#1#0.499918464,#47/150,#entrada de infectados factores externos
    im2=0.001628058, #Probabilidad de infeccion de migrantes despues del 2015
    of=1.000002074,#53, #offset
    rho = 1/(9*52), # influenza tipo A 1/(6*52)  - tipo B (1/(12*365))
    phi = 1/(60*52), # tasa de paso a poblacion mayores de 60 anos
    mud = 0.0573578, # porcentaje de aumento en la tasa de muerte por la influenza
    va=1, #con caida en vacunacion (binaria)
    dist=1,#con caida en vacunacion (binaria)
    rvac=rvac, #caida en cobertura en porcentaje
    tdis=tdis #tiempo de distanciamiento en anos
  )
  y=c(
    S1=44254975*(1-0.00005)-12000, # condiciones iniciales 2009 
    V1P=0,
    V1G=0,
    E1=0,
    I1=100,
    R1=6000,
    S2=44254975*0, # condiciones iniciales 2009 - Proporcion de adultos mayores de 60 anos de acuerdo al censo realizado por el DANE
    V2=0,
    E2=0,
    I2=100,
    R2=6000,
    W=0,
    W1=2028,
    C=0,
    CV1=0,
    CV2=0,
    CV3=0)
  final <-1092#544 # tiempo final de simulacion 60 anos (2030-1066)
  time<-1:final # tiempo para simulacion
  t <- proc.time() # Inicia el cron?metro
  OP <- ode(y,time,flu3,parms,method = "ode45") 
  proc.time()-t
  BASS<-data.frame(OP) 
  
  BASE_Y<-data.frame(semana=1,casos_ci=0,casos_t=0)
  for (i in 2:final) { 
    #Casos
    tpop = OP[i,2]+OP[i,3]+OP[i,4]+OP[i,5]+OP[i,6]
    cas = (BASS$C[i]-BASS$C[i-1])*100000/tpop
    cas0 = (BASS$C[i]-BASS$C[i-1])
    ## Totales
    rec = data.frame(semana=i,casos_ci=cas,casos_t=cas0)
    BASE_Y<-rbind(BASE_Y,rec)
  }
  cas_Y<-NULL
  for (j in 1:21) {
    #j=12
    s1<-BASE_Y$casos_t[(j*52-51):(j*52)]
    s2<-sum(s1)
    re<-data.frame(year=2008+j,cases=s2)
    cas_Y<-rbind(cas_Y,re)
  }
  return = list(c(BASE_Y$casos_t[572:final]),c(BASE_Y$casos_ci[572:final]),BASS[572:final,],cas_Y)# Genera la base de casos del 2019 al 2030 (BASE_Y casos y casos por cien mil) y toda la tabla de los estados del modelo (BASS)
}

######################################
#             Escenarios             #
######################################

# Sin COVID-19 (baseline)
S0<-simuflu(0,0)

# Distanciamiento 6 meses
S6_10<-simuflu(0.5,0.1) # distanciamiento seis meses con caida en cobertura del 10%
S6_25<-simuflu(0.5,0.25) # distanciamiento seis meses con caida en cobertura del 25%
S6_50<-simuflu(0.5,0.5) # distanciamiento seis meses con caida en cobertura del 50%

# Distanciamiento un ano
S1_10<-simuflu(1,0.1) # distanciamiento un ano con caida en cobertura del 10%
S1_25<-simuflu(1,0.25) # distanciamiento un ano con caida en cobertura del 25%
S1_50<-simuflu(1,0.5) # distanciamiento un ano con caida en cobertura del 50%


###########################  Figura casos todo el periodo

par(mfrow=c(1,2))

t=2019:2029
plot(t,S0[[4]][["cases"]][11:21],type="l",lty=1,col="black",lwd=2,xlab="",ylab="Annual cases",main="Influenza cases (2019-2029)",ylim=c(0,max(S0[[4]][["cases"]][10:20])),las=1)
lines(t,S6_10[[4]][["cases"]][11:21],type="l",col="blue",lwd=2)
lines(t,S6_25[[4]][["cases"]][11:21],type="l",col="green",lwd=2)
lines(t,S6_50[[4]][["cases"]][11:21],type="l",col="red",lwd=2)
lines(t,S1_10[[4]][["cases"]][11:21],type="l",col="blue",lwd=2,lty=2)
lines(t,S1_25[[4]][["cases"]][11:21],type="l",col="green",lwd=2,lty=2)
lines(t,S1_50[[4]][["cases"]][11:21],type="l",col="red",lwd=2,lty=2)
legend("topright", legend=c("Baseline","Six months (coverage -10%)  ","Six months (coverage -25%)  ","Six months (coverage -50%)  ",
                            "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.6)

########################### Figura ampliada 2020 al 2022
t2<-2020:2023
plot(t2,S0[[4]][["cases"]][12:15],type="l",lty=1,col="black",lwd=2,xlab="",ylab="",main="Influenza cases (2020-2023)",las=1,ylim=c(min(S1_50[[4]][["cases"]][12:15]),max(S6_50[[4]][["cases"]][12:15])),xaxt = 'n') # make sure to specify yaxt = 'n'
axis(side = 1, at = c(2020,2021,2022,2023,2024), labels = c('2020','2021','2022','2023','2024'))
lines(t2,S6_10[[4]][["cases"]][12:15],type="l",col="blue",lwd=2)
lines(t2,S6_25[[4]][["cases"]][12:15],type="l",col="green",lwd=2)
lines(t2,S6_50[[4]][["cases"]][12:15],type="l",col="red",lwd=2)
lines(t2,S1_10[[4]][["cases"]][12:15],type="l",col="blue",lwd=2,lty=2)
lines(t2,S1_25[[4]][["cases"]][12:15],type="l",col="green",lwd=2,lty=2)
lines(t2,S1_50[[4]][["cases"]][12:15],type="l",col="red",lwd=2,lty=2)


################################################################################
#                               Con proteccion                                 #
################################################################################

# Total personas en el modelo
Total<- S0[[3]][["S1"]]+S0[[3]][["V1P"]]+S0[[3]][["V1G"]]+S0[[3]][["E1"]]+S0[[3]][["I1"]]+S0[[3]][["R1"]]+S0[[3]][["S2"]]+S0[[3]][["V1P"]]+S0[[3]][["V2"]]+S0[[3]][["E2"]]+S0[[3]][["I2"]]+S0[[3]][["R2"]]

##### Proporcion de protegidos sin COVID-19
P01<-(S0[[3]][["V1P"]]+S0[[3]][["V1G"]]+S0[[3]][["V2"]])/Total # INFLUENZA

##### Proporcion de protegidos con 6 meses de distanciamiento con caida en cobertura del 10%
P6_101<-(S6_10[[3]][["V1P"]]+S6_10[[3]][["V1G"]]+S6_10[[3]][["V2"]])/Total # INFLUENZA
##### Proporcion de protegidos con 6 meses de distanciamiento con caida en cobertura del 25%
P6_251<-(S6_25[[3]][["V1P"]]+S6_25[[3]][["V1G"]]+S6_25[[3]][["V2"]])/Total
##### Proporcion de protegidos con 6 meses de distanciamiento con caida en cobertura del 50%
P6_501<-(S6_50[[3]][["V1P"]]+S6_50[[3]][["V1G"]]+S6_50[[3]][["V2"]])/Total

##### Proporcion de protegidos con un ano de distanciamiento con caida en cobertura del 10%
P1_101<-(S1_10[[3]][["V1P"]]+S1_10[[3]][["V1G"]]+S1_10[[3]][["V2"]])/Total
##### Proporcion de protegidos con un ano de distanciamiento con caida en cobertura del 25%
P1_251<-(S1_25[[3]][["V1P"]]+S1_25[[3]][["V1G"]]+S1_25[[3]][["V2"]])/Total
##### Proporcion de protegidos con un ano de distanciamiento con caida en cobertura del 50%
P1_501<-(S1_50[[3]][["V1P"]]+S1_50[[3]][["V1G"]]+S1_50[[3]][["V2"]])/Total

### Figura vacunacion
t=1:521
plot(t/52+2020,P01,type="l",lty=1,col="black",lwd=2,ylim=c(min(P1_501),max(P01)),xlab ="",ylab="Proportion")#,main="Protected by influenza vaccine")
lines(t/52+2020,P6_101,type="l",col="blue",lwd=2)
lines(t/52+2020,P6_251,type="l",col="green",lwd=2)
lines(t/52+2020,P6_501,type="l",col="red",lwd=2)
lines(t/52+2020,P1_101,type="l",col="blue",lwd=2,lty=2)
lines(t/52+2020,P1_251,type="l",col="green",lwd=2,lty=2)
lines(t/52+2020,P1_501,type="l",col="red",lwd=2,lty=2)
legend("bottomright", legend=c("Baseline","Six months (coverage -10%)  ","Six months (coverage -25%)  ","Six months (coverage -50%)  ",
                            "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.8)

####################################################################################################
#                                                                                                  #
#                               Datos de salida del modelo                                         #
#                                                                                                  #
####################################################################################################
# Pronostico casos
SC<-data.frame("id-variable"=rep("UE-INF-SC",11),"Unidades temporales"=2019:2029,Valor_proyectado=S0[[4]][["cases"]][11:21])
CV610<-data.frame("id-variable"=rep("UE-INF-CV610",11),"Unidades temporales"=2019:2029,Valor_proyectado=S6_10[[4]][["cases"]][11:21])
CV625<-data.frame("id-variable"=rep("UE-INF-CV625",11),"Unidades temporales"=2019:2029,Valor_proyectado=S6_25[[4]][["cases"]][11:21])
CV650<-data.frame("id-variable"=rep("UE-INF-CV650",11),"Unidades temporales"=2019:2029,Valor_proyectado=S6_50[[4]][["cases"]][11:21])
CV110<-data.frame("id-variable"=rep("UE-INF-CV110",11),"Unidades temporales"=2019:2029,Valor_proyectado=S1_10[[4]][["cases"]][11:21])
CV125<-data.frame("id-variable"=rep("UE-INF-CV125",11),"Unidades temporales"=2019:2029,Valor_proyectado=S1_25[[4]][["cases"]][11:21])
CV150<-data.frame("id-variable"=rep("UE-INF-CV150",11),"Unidades temporales"=2019:2029,Valor_proyectado=S1_50[[4]][["cases"]][11:21])
BASE_PRO<-rbind(SC,CV610,CV625,CV650,CV110,CV125,CV150)

setwd("~/modelo-epidemiologico")
source("0config.R")
library("WriteXLS")
WriteXLS(BASE_PRO, "INFLUENZA_Pronostico.xlsx")

# Pronostico proporcion de protegidos 
V1SC<-data.frame("id-variable"=rep("UE-INF-PV1SC",521),"Unidades temporales"=t/52+2020,Valor_proyectado=P01)
V1610<-data.frame("id-variable"=rep("UE-INF-PV1610",521),"Unidades temporales"=t/52+2020,Valor_proyectado=P6_101)
V1625<-data.frame("id-variable"=rep("UE-INF-PV1625",521),"Unidades temporales"=t/52+2020,Valor_proyectado=P6_251)
V1650<-data.frame("id-variable"=rep("UE-INF-PV1650",521),"Unidades temporales"=t/52+2020,Valor_proyectado=P6_501)
V1110<-data.frame("id-variable"=rep("UE-INF-PV1110",521),"Unidades temporales"=t/52+2020,Valor_proyectado=P1_101)
V1125<-data.frame("id-variable"=rep("UE-INF-PV1125",521),"Unidades temporales"=t/52+2020,Valor_proyectado=P1_251)
V1150<-data.frame("id-variable"=rep("UE-INF-PV1150",521),"Unidades temporales"=t/52+2020,Valor_proyectado=P1_501)

# Base completa de protegidos
BASE_VAC<-rbind(V1SC,V1610,V1625,V1650,V1110,V1125,V1150)
WriteXLS(BASE_VAC, "INFLUENZA_Pronostico_vac.xlsx")

