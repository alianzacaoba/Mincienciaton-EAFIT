library(beepr)
library(deSolve)

##0. Organizando espacio de trabajo ####
current_directory <- getwd()
print(current_directory)
setwd(current_directory)
setwd("~/modelo-epidemiologico")
source("0config.R")


#################################
#      Modelo parotiditis       #
#################################

mumps <- function(time, y, parms) {
  with(as.list(c(y, parms)),{

    vac1<-W>in1 #  Entrada de vacuna 1 en 1980. Cero hasta 1980 y 1 en adelante
    vac2<-W>in2 # Entrada de vacuna 2 en 1998. Cero hasta 1998 y 1 en adelante
    pdis = (W>18330)*(W<18695)#marzo de 2020 al marzo de 2021 # distanciamiento es 1 entre marzo de 2020 y marzo de 2021

    mu3 = vac2*((cob2-exp(-0.4*(W-9489)/365)))*(1-0.5*pdis*va) # tasa de vacunacion de la vacuna 2
    mu2 = ((cob1*(1-exp(-0.1408442 *(W-3284)/365))))*vac1*(1-0.5*pdis*va)
    
    mu1 = (1-mu2*ef1) # tasa de ingreso de susceptibles 
    
    B1 = (B*cos(2*pi/of*(W))+B)*(1-0.4*pdis*dist)# tasa de trasmision de susceptibles totales
    B2 =B1  # tasa de transmision susceptibles vacunados despues de perder inmunidad
    N <- S1+V1+V2+S2+E+I+R # Numero total de personas
    
    ##### Crecimiento poblacional
    
    n <- ((W*-1.238913e-06 + 3.566146e-02))/365 #funcion de crecimiento Leandro
    mu <- ((0.008251563  +  -4.65767e-07*W + -1.9377e-10 *W^2 + 7.939302e-14*W^3 + -1.160313e-17*W^4 + 8.236164e-22*W^5 + -2.837473e-26*W^6 + 3.794165e-31*W^7))/365
    
    # Estados
    
    ### Entrada de casos de 2015 en adelante por migración positiva
    
    la2<-im*(4.16/1000)/365*N*(W>16790)#*(W<18980)
    
    ######## Ecuaciones 
    
    dS1 <- mu1*N*n - B1*S1*I/N - mu*S1 -la # susceptibles 1
    dV1 <- mu2*N*n*ef1 -tau*V1 -mu3*ef2*V1*(1/(365*5)) -mu*V1 # vacunados V1
    dV2 <- mu3*ef2*V1*(1/(365*5)) -d*V2 -mu*V2 # Vacunados V2
    dS2 <- d*V2 +tau*V1 -B2*S2*I/N -mu*S2 # Susceptibles 2
    dE <- B1*S1*I/N + B2*S2*I/N -al*E -mu*E + la #+ la2 # Expuestos
    dI <- al*E -ga*I -mu*I  + la2 # Infectados
    dR <- ga*I -mu*R # Recuperados
    dW = 1  # tiempo de simulacion # temporizador para ser utilizado como variable W
    dC = B1*S1*I/N + B2*S2*I/N + la + la2 # Nuevos casos # la es la llegada de infectados importados o externos
    list(c(dS1,dV1,dV2,dS2,dE,dI,dR,dW,dC))})
}

################################
#      Simulacio de prueba    #
################################

parms <- c(
  B=0.185759218,#
  cob1 = 0.95, # cobertura vacuna 1
  cob2 = 0.85, # cobertura vacuna 2
  ef1 = 0.88, # efectividad vacuna 1
  ef2 = 0.88, # efectividad vacuna 2
  tau = 0.00034, # tasa de perdida de inmunidad vacuna 1
  d = 0.00034/2, # tasa de perdida de inmunidad vacuna 2
  al = 0.05,# tasa de paso de expuestos a infectados
  ga = 0.167, # tasa de paso de infectados a recuperados
  la = 0.472586574,#probabilidad de entrada de infectados
  of=365, # periodo entre picos de infeccion
  in1 = 365*10, # entrada vacuna 1 1980
  in2 = 365*28, #entra de vacuna 2 1998
  im= 0.001051628,#0.005395314,#0.001131716
  va=0,
  dist=0
)
y=c(
    S1=21.5e6-100, # condiciones iniciales 
    V1=0,
    V2=0,
    S2=0,
    E=0,
    I=1,
    R=0,
    W=0,
    C=0)
final <-365*70 # tiempo final de simulacion 60 anos (2030)
time<-1:final # tiempo para simulacion
t <- proc.time() # Inicia el cronometro
OP <- ode(y,time,mumps,parms,method = "ode2") 
proc.time()-t
BASS<-data.frame(OP) # Base de datos resultados data frame

### Total habitantes
TH=OP[,2]+OP[,3]+OP[,4]+OP[,5]+OP[,6]+OP[,7]+OP[,8]

# Nuevos casos

BASE_Y<-NULL # calcula base por ans de casos por 100.000 habitantes y casos totales
for (i in 1:60) { 
  fin = i*365
  tpop = OP[fin,2]+OP[fin,3]+OP[fin,4]+OP[fin,5]+OP[fin,6]+OP[fin,7]+OP[fin,8]
  cas = (BASS$C[i*365]-BASS$C[i*365-364])*100000/tpop
  cas0 = (BASS$C[i*365]-BASS$C[i*365-364])
  rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0)
  BASE_Y<-rbind(BASE_Y,rec)
}

#########################################################################################
#                                                                                       #
#                                 Calibracion 3 parametros                              #
#                                                                                       #
#########################################################################################

# Datoss de entrada
library('readxl')
setwd(parotiditis_trusted)
casos<- read_excel("mumps.xlsx") # Casos reportados anualmente en la OMS con estimaciï¿½n de casos de la PAHO para 1991 a 1994

## Funcion simular modelo
sirvforc = function(optbeta,la2,im) {
  
  parms <- c(
    B=optbeta,# 
    cob1 = 0.95, # cobertura vacuna 1
    cob2 = 0.85, # cobertura vacuna 2
    ef1 = 0.88, # efectividad vacuna 1
    ef2 = 0.88, # efectividad vacuna 2
    tau = 0.00034, # tasa de pedida de inmunidad vacuna 1
    d = 0.00034/2, # tasa de perdida de inmunidad vacuna 2
    al = 0.05,# tasa de paso de expuestos a infectados
    ga = 0.167, # tasa de paso de infectados a recuperados
    la = la2,#
    of=365, # 
    in1 = 365*10, # entrada vacuna 1 1980
    in2 = 365*28, #entra de vacuna 2 1998
    im = im,
    va=0,
    dist=0
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
  final <-365*50 # tiempo final de simulacion 60 anos (2030)
  time<-1:final # tiempo para simulacion
  t <- proc.time() # Inicia el cronometro
  OP <- ode(y,time,mumps,parms,method = "ode2") 
  proc.time()-t
  BASS<-data.frame(OP) # Base de datos resultados data frame
  
  BASE_Y<-NULL # cicula base por anos de casos por 100.000 habitantes y casos totales
  for (i in 1:50) { 
    fin = i*365
    tpop = OP[fin,2]+OP[fin,3]+OP[fin,4]+OP[fin,5]+OP[fin,6]+OP[fin,7]+OP[fin,8]
    cas = (BASS$C[i*365]-BASS$C[i*365-364])*100000/tpop
    cas0 = (BASS$C[i*365]-BASS$C[i*365-364])
    rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0)
    BASE_Y<-rbind(BASE_Y,rec)
  }
  return = c(BASE_Y$casos_t[22:25],BASE_Y$casos_t[32:50])
}

##### Aqui se hacen las funciones a optimizar

LSmodelo = function(vect1) { #least squares
  loss = sqrt(sum((sirvforc(vect1[1],vect1[2],vect1[3])-casos$Casos)^2))
  return = loss
}

MLENormmodelo = function(vect1) { # Maximun likelihood estimation
  loss = -(sum(dnorm(casos$Casos,mean=sirvforc(vect1[1],vect1[2],vect1[3]),sd=100,log=TRUE)))
  return = loss
}

MLEPoissmodelo = function(vect1) {
  loss = -(sum(dpois(casos$Casos,lambda = sirvforc(vect1[1],vect1[2],vect1[3]) ,log=TRUE)))
  return = loss
}

### A continuacion tenemos las calibraciones para tres parÃ metros

R1<-optim(par=c(1,1,1),fn=LSmodelo,method = "L-BFGS-B",lower = c(0.1,0.1,0), upper = c(0.5,1,0.01))
R2<-optim(par=c(1,1,1), MLENormmodelo,  method = "L-BFGS-B",lower = c(0.1,0.1,0), upper = c(0.5,1,0.01))
#R3<-optim(par=c(1,1,1), fn=MLEPoissmodelo,  method = "L-BFGS-B",lower = c(0.1,0.1,0), upper = c(0.5,1,0.01))

######################################################################
#               Comparison real vs fitting                           #
######################################################################
## From 2001 to 2019
# real
cases_d<-data.frame(year=casos$Ano[5:23],Cases=casos$Casos[5:23],Status=rep("Data",19))
cases_d2<-data.frame(year=casos$Ano[1:23],Cases=casos$Casos[1:23],Status=rep("Data",23))
parms <- c(
  B=0.184490870,# 
  cob1 = 0.95, # 
  cob2 = 0.85, # cobertura vacuna 2
  ef1 = 0.88, # efectividad vacuna 1
  ef2 = 0.88, # efectividad vacuna 2
  tau = 0.00034, # tasa de perdida de inmunidad vacuna 1
  d = 0.00034/2, # tasa de perdida de inmunidad vacuna 2
  al = 0.05,# tasa de paso de expuestos a infectados
  ga = 0.167, # tasa de paso de infectados a recuperados
  la = 0.472586574,#entrada de infectados externos
  of=365, # periodo entre picos de infeccion
  in1 = 365*10, # entrada vacuna 1 1980
  in2 = 365*28, #entra de vacuna 2 1998
  im = 0.001051628,
  va=0,
  dist=0
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
final <-365*50 # tiempo final de simulaciïon 60 anos (2030)
time<-1:final # tiempo para simulacion
t <- proc.time() # Inicia el cronometro
OP <- ode(y,time,mumps,parms,method = "ode2") 
proc.time()-t
BASS<-data.frame(OP) # Base de datos resultados data frame


BASE_Y<-NULL # calcula base por anos de casos por 100.000 habitantes y casos totales
for (i in 1:50) { 
  fin = i*365
  tpop = OP[fin,2]+OP[fin,3]+OP[fin,4]+OP[fin,5]+OP[fin,6]+OP[fin,7]+OP[fin,8]
  cas = (BASS$C[i*365]-BASS$C[i*365-364])*100000/tpop
  cas0 = (BASS$C[i*365]-BASS$C[i*365-364])
  rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0)
  BASE_Y<-rbind(BASE_Y,rec)
}
predic_d<-data.frame(year=BASE_Y$ano[32:50],Cases=BASE_Y$casos_t[32:50],"Status"=rep("Model",19))
BASE_rvsp<-rbind(cases_d,predic_d)
library(ggplot2)
library(dplyr)
ggplot(BASE_rvsp, aes(x=year, y=Cases, group=Status)) +
  geom_line(aes(linetype=Status, color=Status),size=1.5)+ggtitle("Mumps model vs reported cases")

predic_d2<-data.frame(year=BASE_Y$ano[10:50],Cases=BASE_Y$casos_t[10:50],"Status"=rep("Model",41))
BASE_rvsp2<-rbind(cases_d2,predic_d2)
legend_title<-""
ggplot(BASE_rvsp2, aes(x=year, y=Cases, group=Status)) + 
  geom_line(aes(linetype=Status, color=Status, size=Status))+
  geom_point(aes(color=Status, size=Status))+ 
  scale_linetype_manual(values=c("blank","solid"))+
  scale_color_manual(values=c('#F8766D','#00BCD8'))+
  scale_size_manual(values=c(2, 1.8)) + theme(legend.title=element_blank())

######################################################################################
#                                                                                    #
#                         Salida de datos modelo                                     #
#                                                                                    #
######################################################################################
datos<-c(rep(NA,21),casos$Casos[1:4],rep(NA,6),casos$Casos[5:23])
BASE_OUT<-data.frame(Unidades_temporales=BASE_Y$ano,id_variable=rep("UE-MUM-C",50),Valor_proyectado=BASE_Y$casos_t,
                     "Proyectado/100.000 hab"=BASE_Y$casos_ci,"Valor observado"=datos,Error=abs(BASE_Y$casos_t-datos)/datos)
library("WriteXLS")
setwd("~/modelo-epidemiologico")
source("0config.R")
setwd(parotiditis_refined)
WriteXLS(BASE_OUT, "mumps_Observado-Pronostico.xlsx")

