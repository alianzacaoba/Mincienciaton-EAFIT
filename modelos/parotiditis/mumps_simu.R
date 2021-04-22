########################################################################################
#                                                                                      #
#            Modelo de parotiditis calibrado con escenarios de simulacion              #
#                                                                                      #
########################################################################################

library(deSolve)

###############################
#    Modelo compartimentado   #
###############################

mumps4 <- function(time, y, parms) {
  with(as.list(c(y, parms)),{
    
    vac1<-W>in1 #  Entrada de vacuna 1 en 1980. Cero hasta 1980 y 1 en adelante
    vac2<-W>in2 # Entrada de vacuna 2 en 1998. Cero hasta 1998 y 1 en adelante
    pdis = (W>18330)*(W<(18330+tdis*365))#marzo de 2020 al 
    pvac = (W>18330)*(W<(18330+tdis*365))#marzo de 2020 al 
    
    mu3 = vac2*((cob2-exp(-0.4*(W-9489)/365)))*(1-rvac*pdis*va) # tasa de vacunacion de la vacuna 2
    mu2 = ((cob1*(1-exp(-0.1408442 *(W-3284)/365))))*vac1*(1-rvac*pvac*va)
    
    mu1 = (1-mu2*ef1) # tasa de ingreso de susceptibles 
    
    B1 = (B*cos(2*pi/of*(W))+B)*(1-0.4*pdis*dist)# tasa de trasmision de susceptibles totales
    B2 =B1  # tasa de transmision susceptibles vacunados despues de perder inmunidad
    N <- S1+V1+V2+S2+E+I+R # Numero total de personas
    
    ##### Crecimiento poblacional
    
    n <- ((W*-1.238913e-06 + 3.566146e-02))/365 #funcion de crecimiento 
    mu <- ((0.008251563  +  -4.65767e-07*W + -1.9377e-10 *W^2 + 7.939302e-14*W^3 + -1.160313e-17*W^4 + 8.236164e-22*W^5 + -2.837473e-26*W^6 + 3.794165e-31*W^7))/365
    
    # Estados
    
    ### Entrada de casos de 2015 en adelante (solo para parotiditis)
    
    la2<-im*(4.16/1000)/365*N*(W>16790)*(W<18980)
    
    ######## Ecuaciones 
    
    dS1 <- mu1*N*n - B1*S1*I/N - mu*S1 -la # susceptibles 1
    dV1 <- mu2*N*n*ef1 -tau*V1 -mu3*ef2*V1*(1/(365*5)) -mu*V1 # vacunados V1
    dV2 <- mu3*ef2*V1*(1/(365*5)) -d*V2 -mu*V2 # Vacunados V2
    dS2 <- d*V2 +tau*V1 -B2*S2*I/N -mu*S2 # Susceptibles 2
    dE <- B1*S1*I/N + B2*S2*I/N -al*E -mu*E + la*(1-pdis) #+ la2 # Expuestos
    dI <- al*E -ga*I -mu*I  + la2 # Infectados
    dR <- ga*I -mu*R # Recuperados
    
    ############## Adicionales para calculos
    
    dW = 1  # tiempo de simulacion # temporizador para ser utilizado como variable W
    dC = B1*S1*I/N + B2*S2*I/N + la + la2 # Nuevos casos # la es la llegada de infectados importados o externos
    
    list(c(dS1,dV1,dV2,dS2,dE,dI,dR,dW,dC))})
}

######################################
## Funcion para simular escenarios   #
######################################
# Pide tiempo de distanciamiento (tdis)
# y porcentaje de caida en la cobertura
# de vacunacion durante la pandemia (rvac)

simumumps = function(tdis,rvac) {
  
  parms <- c(
    B=0.184490870,#
    cob1 = 0.95, # cobertura vacuna 1
    cob2 = 0.85, # cobertura vacuna 2
    ef1 = 0.88, # efectividad vacuna 1
    ef2 = 0.88, # efectividad vacuna 2
    tau = 0.00034, # tasa de perdida de inmunidad vacuna 1
    d = 0.00034/2, # tasa de perdida de inmunidad vacuna 2
    al = 0.05,# tasa de paso de expuestos a infectados
    ga = 0.167, # tasa de paso de infectados a recuperados
    la = 0.472586574,#probabilidad de entrada de infectados externos
    of=365, # periódo entre picos de infección
    in1 = 365*10, # entrada vacuna 1 1980
    in2 = 365*28, #entra de vacuna 2 1998
    im =  0.001051628,
    va=1,
    dist=1,
    tdis=tdis,
    rvac=rvac
  )
  y=c(
    S1=21.5e6-100, # condiciones iniciales 
    V1=0,
    V2=0,
    S2=0,
    E=0,
    I=100,
    R=0,
    W=0,
    C=0)
  final <-365*60 # tiempo final de simulación 60 anos (2030)
  time<-1:final # tiempo para simulacion
  OP <- ode(y,time,mumps4,parms,method = "ode2") 
  BASS<-data.frame(OP) # Base de datos resultados data frame
  
  
  BASE_Y<-NULL # calcula base por años de casos por 100.000 habitantes y casos totales
  for (i in 1:60) { 
    fin = i*365
    #Casos
    tpop = OP[fin,2]+OP[fin,3]+OP[fin,4]+OP[fin,5]+OP[fin,6]+OP[fin,7]+OP[fin,8]
    cas = (BASS$C[i*365]-BASS$C[i*365-364])*100000/tpop
    cas0 = (BASS$C[i*365]-BASS$C[i*365-364])
    ## Totales
    rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0)
    BASE_Y<-rbind(BASE_Y,rec)
  }
  return = list(c(BASE_Y$casos_t[50:60]),c(BASE_Y$casos_ci[50:60]),BASS)# Genera la base de casos del 2019 al 2020 (BASE_Y casos y casos por cien mil) y toda la tabla de los estados del modelo (BASS)
}

######################################
#             Escenarios             #
######################################

# Sin COVID-19 (baseline)
S0<-simumumps(0,0)

# Distanciamiento 6 meses
S6_10<-simumumps(0.5,0.1) # distanciamiento seis meses con caida en cobertura del 10%
S6_25<-simumumps(0.5,0.25) # distanciamiento seis meses con caida en cobertura del 25%
S6_50<-simumumps(0.5,0.5) # distanciamiento seis meses con caida en cobertura del 50%

# Distanciamiento un año
S1_10<-simumumps(1,0.1) # distanciamiento un ano con caida en cobertura del 10%
S1_25<-simumumps(1,0.25) # distanciamiento un ano con caida en cobertura del 25%
S1_50<-simumumps(1,0.5) # distanciamiento un ano con caida en cobertura del 50%

###########################  Figura casos todo el periodo

par(mfrow=c(1,2))

plot(2019:2029,S0[[1]],type="l",lty=1,col="black",lwd=2,xlab="",ylab="Annual cases",main="Mumps cases (2019-2029)",ylim=c(0,max(S0[[1]])),las=1,xaxt="none")
axis(1, seq(2019,2029,1))
lines(2019:2029,S6_10[[1]],type="l",col="blue",lwd=2)
lines(2019:2029,S6_25[[1]],type="l",col="green",lwd=2)
lines(2019:2029,S6_50[[1]],type="l",col="red",lwd=2)
lines(2019:2029,S1_10[[1]],type="l",col="blue",lwd=2,lty=2)
lines(2019:2029,S1_25[[1]],type="l",col="green",lwd=2,lty=2)
lines(2019:2029,S1_50[[1]],type="l",col="red",lwd=2,lty=2)
legend("topright", legend=c("Baseline","Six months (coverage -10%)  ","Six months (coverage -25%)  ","Six months (coverage -50%)  ",
                            "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.8)

########################### Figura ampliada 2024 en adelante de casos por 100.000 habitantes

plot(2024:2029,S0[[1]][6:11],type="l",lty=1,col="black",lwd=2,xlab="",ylab="Annual cases",main="Mumps cases (2024-2029)",ylim=c(min(S0[[1]][6:11]),max(S1_50[[1]][6:11])),las=1)
lines(2024:2029,S6_10[[1]][6:11],type="l",col="blue",lwd=2)
lines(2024:2029,S6_25[[1]][6:11],type="l",col="green",lwd=2)
lines(2024:2029,S6_50[[1]][6:11],type="l",col="red",lwd=2)
lines(2024:2029,S1_10[[1]][6:11],type="l",col="blue",lwd=2,lty=2)
lines(2024:2029,S1_25[[1]][6:11],type="l",col="green",lwd=2,lty=2)
lines(2024:2029,S1_50[[1]][6:11],type="l",col="red",lwd=2,lty=2)

################################################################################
#                               Con proteccion                                 #
################################################################################

# Total personas en el modelo
Total<- (S0[[3]][["S1"]]+S0[[3]][["V1"]]+S0[[3]][["V2"]]+S0[[3]][["S2"]]+S0[[3]][["E"]]+S0[[3]][["I"]]+S0[[3]][["R"]])

##### Proporcion de protegidos sin COVID-19
P01<-S0[[3]][["V1"]]/Total # MMR1
P02<-S0[[3]][["V2"]]/Total # MMR2

##### Proporcion de protegidos con 6 meses de distanciamiento con caida en cobertura del 10%
P6_101<-S6_10[[3]][["V1"]]/Total # MMR1
P6_102<-S6_10[[3]][["V2"]]/Total # MMR2
##### Proporcion de protegidos con 6 meses de distanciamiento con caida en cobertura del 25%
P6_251<-S6_25[[3]][["V1"]]/Total # MMR1
P6_252<-S6_25[[3]][["V2"]]/Total # MMR2
##### Proporcion de protegidos con 6 meses de distanciamiento con caida en cobertura del 50%
P6_501<-S6_50[[3]][["V1"]]/Total # MMR1
P6_502<-S6_50[[3]][["V2"]]/Total # MMR2

##### Proporcion de protegidos con un año de distanciamiento con caida en cobertura del 10%
P1_101<-S1_10[[3]][["V1"]]/Total # MMR1
P1_102<-S1_10[[3]][["V2"]]/Total # MMR2
##### Proporcion de protegidos con un año de distanciamiento con caida en cobertura del 25%
P1_251<-S1_25[[3]][["V1"]]/Total # MMR1
P1_252<-S1_25[[3]][["V2"]]/Total # MMR2
##### Proporcion de protegidos con un año de distanciamiento con caida en cobertura del 50%
P1_501<-S1_50[[3]][["V1"]]/Total # MMR1
P1_502<-S1_50[[3]][["V2"]]/Total # MMR2

#### Tiempo en dias entre el 2019 y 2029
Te=18250:21900

##################################### Figura proporcion de protegidos

par(mfrow=c(1,2))

plot(Te/365+1969,P01[18250:21900],type="l",lty=1,col="black",lwd=2,ylim=c(min(P1_501[18250:21900]),max(P01[18250:21900])),xlab ="",ylab="Proportion",main="Protected by MMR 1 (mumps)")
lines(Te/365+1969,P6_101[18250:21900],type="l",col="blue",lwd=2)
lines(Te/365+1969,P6_251[18250:21900],type="l",col="green",lwd=2)
lines(Te/365+1969,P6_501[18250:21900],type="l",col="red",lwd=2)
lines(Te/365+1969,P1_101[18250:21900],type="l",col="blue",lwd=2,lty=2)
lines(Te/365+1969,P1_251[18250:21900],type="l",col="green",lwd=2,lty=2)
lines(Te/365+1969,P1_501[18250:21900],type="l",col="red",lwd=2,lty=2)
legend("topright", legend=c("Baseline","Six months (coverage -10%)  ","Six months (coverage -25%)  ","Six months (coverage -50%)  ",
                               "One year (coverage -10%)","One year (coverage -25%)","One year (coverage -50%)"),
       col=c("black","blue","green","red","blue","green","red"), lty=c(1,1,1,1,2,2,2),lwd=c(2,2,2,2,2,2,2), cex=0.8)


plot(Te/365+1969,P02[18250:21900],type="l",lty=1,col="black",lwd=2,ylim=c(min(P1_502[18250:21900]),max(P02[18250:21900])),xlab ="",ylab="Proportion",main="Protected by MMR 2 (mumps)")
lines(Te/365+1969,P6_102[18250:21900],type="l",col="blue",lwd=2)
lines(Te/365+1969,P6_252[18250:21900],type="l",col="green",lwd=2)
lines(Te/365+1969,P6_502[18250:21900],type="l",col="red",lwd=2)
lines(Te/365+1969,P1_102[18250:21900],type="l",col="blue",lwd=2,lty=2)
lines(Te/365+1969,P1_252[18250:21900],type="l",col="green",lwd=2,lty=2)
lines(Te/365+1969,P1_502[18250:21900],type="l",col="red",lwd=2,lty=2)

####################################################################################################
#                                                                                                  #
#                               Datos de salida del modelo                                         #
#                                                                                                  #
####################################################################################################
# Pronostico casos
SC<-data.frame("id-variable"=rep("UE-MUM-SC",11),"Unidades temporales"=2019:2029,Valor_proyectado=S0[[1]])
CV610<-data.frame("id-variable"=rep("UE-MUM-CV610",11),"Unidades temporales"=2019:2029,Valor_proyectado=S6_10[[1]])
CV625<-data.frame("id-variable"=rep("UE-MUM-CV625",11),"Unidades temporales"=2019:2029,Valor_proyectado=S6_25[[1]])
CV650<-data.frame("id-variable"=rep("UE-MUM-CV650",11),"Unidades temporales"=2019:2029,Valor_proyectado=S6_50[[1]])
CV110<-data.frame("id-variable"=rep("UE-MUM-CV110",11),"Unidades temporales"=2019:2029,Valor_proyectado=S1_10[[1]])
CV125<-data.frame("id-variable"=rep("UE-MUM-CV125",11),"Unidades temporales"=2019:2029,Valor_proyectado=S1_25[[1]])
CV150<-data.frame("id-variable"=rep("UE-MUM-CV150",11),"Unidades temporales"=2019:2029,Valor_proyectado=S1_50[[1]])
BASE_PRO<-rbind(SC,CV610,CV625,CV650,CV110,CV125,CV150)

setwd("~/modelo-epidemiologico")
source("0config.R")
setwd(parotiditis_refined)
library("WriteXLS")
WriteXLS(BASE_PRO, "mumps_Pronostico.xlsx")

# Pronostico proporcion de protegidos para la parotiditis (vacuna MMR1) 
V1SC<-data.frame("id-variable"=rep("UE-MUM-PV1SC",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P01[18250:21900])
V1610<-data.frame("id-variable"=rep("UE-MUM-PV1610",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P6_101[18250:21900])
V1625<-data.frame("id-variable"=rep("UE-MUM-PV1625",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P6_251[18250:21900])
V1650<-data.frame("id-variable"=rep("UE-MUM-PV1650",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P6_501[18250:21900])
V1110<-data.frame("id-variable"=rep("UE-MUM-PV1110",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P1_101[18250:21900])
V1125<-data.frame("id-variable"=rep("UE-MUM-PV1125",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P1_251[18250:21900])
V1150<-data.frame("id-variable"=rep("UE-MUM-PV1150",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P1_501[18250:21900])
# Pronostico proporcion de protegidos para la parotiditis (vacuna MMR2) 
V2SC<-data.frame("id-variable"=rep("UE-MUM-PV2SC",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P02[18250:21900])
V2610<-data.frame("id-variable"=rep("UE-MUM-PV2610",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P6_102[18250:21900])
V2625<-data.frame("id-variable"=rep("UE-MUM-PV2625",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P6_252[18250:21900])
V2650<-data.frame("id-variable"=rep("UE-MUM-PV2650",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P6_502[18250:21900])
V2110<-data.frame("id-variable"=rep("UE-MUM-PV2110",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P1_102[18250:21900])
V2125<-data.frame("id-variable"=rep("UE-MUM-PV2125",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P1_252[18250:21900])
V2150<-data.frame("id-variable"=rep("UE-MUM-PV2150",3651),"Unidades temporales"=Te/365+1969,Valor_proyectado=P1_502[18250:21900])
# Base completa de protegidos
BASE_VAC<-rbind(V1SC,V1610,V1625,V1650,V1110,V1125,V1150,V2SC,V2610,V2625,V2650,V2110,V2125,V2150)
WriteXLS(BASE_VAC, "mumps_Pronostico_vac.xlsx")
