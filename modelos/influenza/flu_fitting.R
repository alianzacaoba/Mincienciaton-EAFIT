#############################################################################
#                                                                           #
#                             Influenza con tres grupos                     #
#                                                                           #
#############################################################################
library(deSolve)
flu2 <- function(time, y, parms) {
  with(as.list(c(y, parms)),{
    
    ###############
    pdis = (W>585)*(W<(585+tdis*52))#marzo de 2020 al 
    pvac = (W>585)*(W<(585+tdis*52))#marzo de 2020 al 
    
    W2=W
    mu2 = ((4.347058e-01 -5.035106e-03*W2 + 2.982051e-05*W2^2 + 7.823542e-07*W2^3 -8.044537e-09*W2^4 + 2.997302e-11*W2^5 -5.035786e-14*W2^6 + 3.212231e-17*W2^7)*(W2<(9*52)) + 0.7058571*(W2>(9*52-1)))*(1-rvac*pvac*va) # Función de cobertura del 2009 en adelante
    mu3 = (0.6166667*(W2>(52*4)))*(1-rvac*pdis*va) # Cobertura de vacunación de grupo de riesgo 2 (gestantes) del 2009 en adelante
    mu4 = 0.5*(1-rvac*pvac*va)#((2.225411e-01 +  3.111241e-02*W2 -6.578987e-04*W2^2 + 4.361894e-06*W2^3 -9.084438e-09*W2^4 -3.937422e-12*W2^5 + 2.165713e-14*W2^6)*(W2<(7*52)) + 0.71*(W2>(7*52-1)))*(1-rvac*pdis*va)# cobertura de vacunacion grupo de riesgo 3 (adultos mayores)
    
    B1 = (B*cos(2*pi/52*(W-of))+B)*(1-0.4*pdis*dist)# tasa de trasmision de susceptibles totales
    N <- S1+V1P+V1G+E1+I1+R1+S2+V2+E2+I2+R2 # Número total de personas
    
    ############# Entrada de casos externos migracion después del 2015 (0 migrantes por 100 mil) que es positiva hasta llegar a 4.16 migrantes por 100 mil en 2020
    la2<-im2*((W*0.016-4.992)/1000)/52*N*(W>312)*(W<624)
    
    ##### Crecimiento poblacional
    
    n <- ((W1*-1.238913e-06 + 3.566146e-02))/52 #función de crecimiento poblacional
    mu <-((0.008251563  +  -4.65767e-07*W1 + -1.9377e-10 *W1^2 + 7.939302e-14*W1^3 + -1.160313e-17*W1^4 + 8.236164e-22*W1^5 + -2.837473e-26*W1^6 + 3.794165e-31*W1^7))/52 #función de tasa de muerte
    
    # Estados
    
    dS1 <- N*n*(1-mu2*ef1) +rho*R1 +tau*V1P +tau*V1G -B1*S1*(I1+I2)/N - N*n*mu3*ef2 - phi*S1 - mu*S1 - im # susceptibles antes de los 60
    dV1P <- N*n*mu2*ef1 -tau*V1P -mu*V1P # vacunados grupo 1 (niños menores a un ano)
    dV1G <- N*n*mu3*ef2 -tau*V1G -mu*V1G # vacunados grupo 2 (madres gestantes)
    dE1 <-  B1*S1*(I1+I2)/N - al*E1 -mu*E1 +im -phi*E1# Expuestos menores de 60 anos
    dI1 <-  al*E1 -ga*I1 -mu*I1 -mud*I1 +la2*(1-0.066) -phi*I1# Infectados menores de 60 anos
    dR1 <-  ga*I1 -rho*R1 -mu*R1 -phi*R1# Recuperados menores de 60 anos
    dS2 <- phi*S1 +rho*R2 +tau*V2 -B1*S2*(I1+I2)/N - mu4*ef3*S2 - mu*S2 -im# susceptibles mayores de 60
    dV2 <- mu4*ef3*S2 -tau*V2 -mu*V2 # vacunados grupo 3 (mayores de 60 anos)
    dE2 <- B1*S2*(I1+I2)/N -al*E2 -mu*E2 +im +phi*E1# Expuestos menores de 60 anos
    dI2 <- al*E2 -ga*I2 -mu*I2 -mud*I2 + la2*0.066 +phi*I1 # Infectados menores de 60 anos
    dR2 <- ga*I2 -rho*R2 -mu*R2 -phi*R1# Recuperados menores de 60 anos
    
    ############## Adicionales para calculos
    
    dW = 1 # inicio en 2009  # tiempo de simulación # temporizador para funciones de vacunacion y de tasa de transmision
    dW1 = 1 # temporizador para funciones de tasa de nacimiento y de mortalidad
    dC = al*(E1+E2) +la2 # Nuevos casos # transmission (al*E) e importados (im)
    dCV1 = N*n*mu2*ef1
    dCV2 = N*n*mu3*ef1
    dCV3 = mu4*ef1*S2
    
    list(c(dS1,dV1P,dV1G,dE1,dI1,dR1,dS2,dV2,dE2,dI2,dR2,dW,dW1,dC,dCV1,dCV2,dCV3))})}
parms <- c(
  B=0.732169990,#
  ef1 = 0.74, # efectividad vacuna pop1
  ef2 = 0.54, # efectividad vacuna pop2
  ef3 = 0.4, # efectividad vacuna pop3
  tau = 1/52, # tasa de pérdida de inmunidad vacuna 
  al = 1/(1/7),# tasa de paso de expuestos a infectados
  ga = 1/(5/7), # tasa de paso de infectados a recuperados
  im = 1.905911899,#1#entrada de infectados factores externos
  im2=0.001628058, #Probabilidad de infección de migrantes despues del 2015
  of=1.000002074,#53, #offset
  rho = 1/(9*52), # influenza tipo A 1/(6*52)  - tipo B (1/(12*365))
  phi = 1/(60*52), # tasa de paso a población mayores de 60 anos
  mud = 0.0573578, # porcentaje de aumento en la tasa de muerte por la influenza
  va=1, #con caida en vacunacion (binaria)
  dist=1,#con caida en vacunacion (binaria)
  rvac=0.5, #caida en cobertura en porcentaje
  tdis=1 #tiempo de distanciamiento en anos
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
final <-572 # (544 hasta el 2019)tiempo final de simulacion 60 anos (2030-1066)
time<-1:final # tiempo para simulacion
t <- proc.time() # Inicia el cronometro
OP <- ode(y,time,flu2,parms,method = "ode45") 
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
for (j in 1:11) {
  s1<-BASE_Y$casos_t[(j*52-51):(j*52)]
  s2<-sum(s1)
  re<-data.frame(year=2008+j,cases=s2)
  cas_Y<-rbind(cas_Y,re)
}
plot(BASS$I1+BASS$I2,type="l")

#############################################
#                                           #
#             Model fitting                 #
#                                           #
#############################################
setwd("~/modelo-epidemiologico")
source("0config.R")
setwd(influenza_trusted)
library("readxl")
casos<- read_excel("FluNetInteractiveReport.xlsx")
new_cases<-data.frame(year=as.numeric(casos$Year),week=as.numeric(casos$Week),positive_cases=as.numeric(casos$`Total number of influenza positive viruses`))
# Casos anuales desde el 2009
ANN<-NULL
for (i in 2009:2019) {
  s1<-subset(new_cases,new_cases$year==i)
  s2<-sum(na.omit(s1$positive_cases))
  rc<-data.frame(year=i,cases=s2)
  ANN<-rbind(ANN,rc)
}
########## Muertes por neumonia
my_data <- read.delim("death2_1998_2008.txt",sep =';')
neumonia<-subset(my_data,my_data$C_DIR1=="J440")
NEU<-NULL
for (i in 1998:2018) {
  sub<-subset(neumonia,neumonia$ANO==i)
  cas<-nrow(sub)
  rec<-data.frame(ANO=i,CASOS=cas)
  NEU<-rbind(NEU,rec)
}
#####
ANN2<-data.frame(ANN,deaths=c(NEU$CASOS[12:21],NA),pro_death=c(NEU$CASOS[12:21],NA)/ANN$cases)
### Tasa de muerte 
mud<-mean(na.omit(ANN2$pro_death))
# Modelo
# flu2

flu3 <- function(time, y, parms) {
  with(as.list(c(y, parms)),{
    
    ###############
    pdis = (W>585)*(W<(585+tdis*52))#marzo de 2020 al 
    pvac = (W>585)*(W<(585+tdis*52))#marzo de 2020 al 
    
    W2=W
    mu2 = ((4.347058e-01 -5.035106e-03*W2 + 2.982051e-05*W2^2 + 7.823542e-07*W2^3 -8.044537e-09*W2^4 + 2.997302e-11*W2^5 -5.035786e-14*W2^6 + 3.212231e-17*W2^7)*(W2<(9*52)) + 0.7058571*(W2>(9*52-1)))*(1-rvac*pvac*va) # Función de cobertura del 2009 en adelante
    mu3 = (0.6166667*(W2>(52*4)))*(1-rvac*pdis*va) # Cobertura de vacunacion de grupo de riesgo 2 (gestantes) del 2009 en adelante
    mu4 = 0.5*(1-rvac*pvac*va)#((2.225411e-01 +  3.111241e-02*W2 -6.578987e-04*W2^2 + 4.361894e-06*W2^3 -9.084438e-09*W2^4 -3.937422e-12*W2^5 + 2.165713e-14*W2^6)*(W2<(7*52)) + 0.71*(W2>(7*52-1)))*(1-rvac*pdis*va)# cobertura de vacunacion grupo de riesgo 3 (adultos mayores)
    
    B1 = (B*cos(2*pi/52*(W-of))+B)*(1-0.4*pdis*dist)# tasa de trasmissión de susceptibles totales
    N <- S1+V1P+V1G+E1+I1+R1+S2+V2+E2+I2+R2 # Numero total de personas
    
    ############# Entrada de casos externos migracion despuée del 2015 (0 migrantes por 100 mil) que es positiva hasta llegar a 4.16 migrantes por 100 mil en 2020
    la2<-im2*((W*0.016-4.992)/1000)/52*N*(W>312)*(W<624)
    
    ##### Crecimiento poblacional
    
    n <- ((W1*-1.238913e-06 + 3.566146e-02))/52 #función de crecimiento poblacional
    mu <-((0.008251563  +  -4.65767e-07*W1 + -1.9377e-10 *W1^2 + 7.939302e-14*W1^3 + -1.160313e-17*W1^4 + 8.236164e-22*W1^5 + -2.837473e-26*W1^6 + 3.794165e-31*W1^7))/52 #funcion de tasa de muerte
    
    # Estados
    
    dS1 <- N*n*(1-mu2*ef1) +rho*R1 +tau*V1P +tau*V1G -B1*S1*(I1+I2)/N - N*n*mu3*ef2 - phi*S1 - mu*S1 - im # susceptibles antes de los 60
    dV1P <- N*n*mu2*ef1 -tau*V1P -mu*V1P # vacunados grupo 1 (niños menores a un ano)
    dV1G <- N*n*mu3*ef2 -tau*V1G -mu*V1G # vacunados grupo 2 (madres gestantes)
    dE1 <-  B1*S1*(I1+I2)/N - al*E1 -mu*E1 +im -phi*E1# Expuestos menores de 60 anos
    dI1 <-  al*E1 -ga*I1 -mu*I1 -mud*I1 +la2*(1-0.066) -phi*I1# Infectados menores de 60 anos
    dR1 <-  ga*I1 -rho*R1 -mu*R1 -phi*R1# Recuperados menores de 60 anos
    dS2 <- phi*S1 +rho*R2 +tau*V2 -B1*S2*(I1+I2)/N - mu4*ef3*S2 - mu*S2 -im# susceptibles mayores de 60
    dV2 <- mu4*ef3*S2 -tau*V2 -mu*V2 # vacunados grupo 3 (mayores de 60 anos)
    dE2 <- B1*S2*(I1+I2)/N -al*E2 -mu*E2 +im +phi*E1# Expuestos menores de 60 anos
    dI2 <- al*E2 -ga*I2 -mu*I2 -mud*I2 + la2*0.066 +phi*I1 # Infectados menores de 60 anos
    dR2 <- ga*I2 -rho*R2 -mu*R2 -phi*R1# Recuperados menores de 60 anos
    
    ############## Adicionales para calculos
    
    dW = 1 # inicio en 2009  # tiempo de simulación # temporizador para funciones de vacunación y de tasa de transmisión
    dW1 = 1 # temporizador para funciones de tasa de nacimiento y de mortalidad
    dC = al*(E1+E2)+la2 # Nuevos casos # transmission (al*E) e importados (im)
    dCV1 = N*n*mu2*ef1
    dCV2 = N*n*mu3*ef1
    dCV3 = mu4*ef1*S2
    
    list(c(dS1,dV1P,dV1G,dE1,dI1,dR1,dS2,dV2,dE2,dI2,dR2,dW,dW1,dC,dCV1,dCV2,dCV3))})}

#Funcion salida de casos
sirvforc = function(optbeta,im,of,im2) {
  parms <- c(
    B=optbeta,#
    ef1 = 0.74, # efectividad vacuna pop1
    ef2 = 0.54, # efectividad vacuna pop2
    ef3 = 0.4, # efectividad vacuna 
    tau = 1/52, # tasa de pérdida de inmunidad vacuna 
    al = 1/(1/7),# tasa de paso de expuestos a infectados
    ga = 1/(5/7), # tasa de paso de infectados a recuperados
    im = im,#1#0.499918464,#47/150,#entrada de infectados factores externos
    im2=im2, #Probabilidad de infección de migrantes después del 2015
    of=of,#53, # periódo entre picos de infección
    rho = 1/(9*53), # influenza tipo A - tipo B 
    phi = 1/(60*52), # tasa de paso a poblacion mayores de 60 anos
    mud = 0.0573578, # porcentaje  de muerte por la influenza
    va=0, #con caida en vacunación (binaria)
    dist=0,#con caida en vacunación (binaria)
    rvac=0, #caída en cobertura en porcentaje
    tdis=0 #tiempo de distanciamiento en anos
  )
  y=c(
    S1=44254975*(1-0.00005)-12000, # condiciones iniciales 2009 
    V1P=0,
    V1G=0,
    E1=0,
    I1=100,
    R1=6000,
    S2=44254975*0, 
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
  final <-572 # (2009-2019)
  time<-1:final # tiempo para simulacion
  t <- proc.time() # Inicia el cronometro
  OP <- ode(y,time,flu2,parms,method = "ode45") 
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
  for (j in 1:11) {
    s1<-BASE_Y$casos_t[(j*52-51):(j*52)]
    s2<-sum(s1)
    re<-data.frame(year=2008+j,cases=s2)
    cas_Y<-rbind(cas_Y,re)
  }
  return = cas_Y$cases
}

#Funcion salida de casos (sin vacunas)
sirvforc2 = function(optbeta,im,of,im2) {
  parms <- c(
    B=optbeta,#
    ef1 = 0.74, # efectividad vacuna pop1
    ef2 = 0.54, # efectividad vacuna pop2
    ef3 = 0.4, # efectividad vacuna pop3
    tau = 1/52, # tasa de perdida de inmunidad vacuna 
    al = 1/(1/7),# tasa de paso de expuestos a infectados
    ga = 1/(5/7), # tasa de paso de infectados a recuperados
    im = im,#1#0.499918464,#47/150,#entrada de infectados factores externos
    im2=im2, #Probabilidad de infección de migrantes después del 2015
    of=of,#53, # periodo entre picos de infeccion
    rho = 1/(9*53), # influenza 
    phi = 1/(60*52), # tasa de paso a poblacion mayores de 60 anos
    mud = 0.0573578, # porcentaje  de muerte por la influenza
    va=0, #con caida en vacunacion (binaria)
    dist=0,#con caida en vacunacion (binaria)
    rvac=0, #caida en cobertura en porcentaje
    tdis=0 #tiempo de distanciamiento en anos
  )
  y=c(
    S1=44254975*(1-0.00005)-12000, # condiciones iniciales 2009 
    V1P=0,
    V1G=0,
    E1=0,
    I1=100,
    R1=6000,
    S2=44254975*0, 
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
  final <-572 # (2009-2019)
  time<-1:final # tiempo para simulacion
  t <- proc.time() # Inicia el cronometro
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
  for (j in 1:11) {
    s1<-BASE_Y$casos_t[(j*52-51):(j*52)]
    s2<-sum(s1)
    re<-data.frame(year=2008+j,cases=s2)
    cas_Y<-rbind(cas_Y,re)
  }
  return = cas_Y$cases
}

##### Aqui se hacen las funciones para poner en optim()
LSmodelo = function(vect1) {#least squares
  loss = sqrt(sum((sirvforc(vect1[1],vect1[2],vect1[3],vect1[4])-ANN$cases)^2))
  return = loss
}
### A continuacion tenemos las calibraciones para tres parÃ metros
R1<-optim(par=c(1,1,1,1),fn=LSmodelo,method = "L-BFGS-B",lower = c(0.7,1/52,1,0.0000001), upper = c(2,2,52,0.01))
#Solve 0.732169990 1.905911899 1.000002074 0.001628058
load(file="sintoflu_2_12.Rdata")
completo<-sirvforc(R1$par[1],R1$par[2],R1$par[3],R1$par[4])
sin_vacunas<-sirvforc2(R1$par[1],R1$par[2],R1$par[3],R1$par[4])
#sa<-sirvforc(1.197352642,2,2,0.002)
plot(2009:2019,completo,type="l",ylim=c(0,max(max(completo),max(ANN$cases))),col="cadetblue2",lwd=3,xlab=" ",ylab="Cases per year")

points(2009:2019,ANN$cases,col="red",lwd=1,pch=16)
legend("topright", legend=c("Model  ","Data  "), pch=c(NA,16),
       col=c("cadetblue2","red"), lty=c(1,0),lwd=c(2,1),cex=0.8)
#############################
parms <- c(
  B=0.732169990,#
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
  rvac=0.5, #caida en cobertura en porcentaje
  tdis=1 #tiempo de distanciamiento en anos
)
y=c(
  S1=44254975*(1-0.00005)-12000, # condiciones iniciales 2009 
  V1P=0,
  V1G=0,
  E1=0,
  I1=100,
  R1=6000,
  S2=44254975*0, 
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
final <-572 # (544 hasta el 2019)tiempo final de simulación 60 anos (2030-1066)
time<-1:final # tiempo para simulacion
t <- proc.time() # Inicia el cronometro
OP <- ode(y,time,flu2,parms,method = "ode45") 
proc.time()-t
BASS<-data.frame(OP) 

BASE_Y<-NULL # calcula base por anos de casos por 100.000 habitantes y casos totales
for (i in 1:11) { 
  fin = i*52
  tpop = OP[fin,2]+OP[fin,3]+OP[fin,4]+OP[fin,5]+OP[fin,6]+OP[fin,7]+OP[fin,8]+OP[fin,9]+OP[fin,10]+OP[fin,11]+OP[fin,12]
  cas = (BASS$C[i*52]-BASS$C[i*52-51])*100000/tpop
  cas0 = (BASS$C[i*52]-BASS$C[i*52-51])
  rec = data.frame(ano=i+2008,casos_ci=cas,casos_t=cas0)
  BASE_Y<-rbind(BASE_Y,rec)
}
######################################################################################
#                                                                                    #
#                         Salida de datos modelo                                     #
#                                                                                    #
######################################################################################
datos<-ANN$cases
BASE_OUT<-data.frame(Unidades_temporales=2009:2019,id_variable=rep("UE-FLU-C",11),Valor_proyectado=BASE_Y$casos_t,
                     "Proyectado/100.000 hab"=BASE_Y$casos_ci,"Valor observado"=datos,Error=abs(BASE_Y$casos_t-datos)/datos)
setwd("~/modelo-epidemiologico")
source("0config.R")
setwd(influenza_refined)
library("WriteXLS")
WriteXLS(BASE_OUT, "flu_Observado-Pronostico.xlsx")

