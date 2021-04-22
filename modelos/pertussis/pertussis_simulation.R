# calibracion MLE 

library(deSolve)
library(dfoptim)
options(scipen=999)

############### Aquí se cargan los datos de contagios ############################
#### Los datos que se suministran aquí son datos de prueba para la calibracion, tomados del correo del dr Ange, y vienen
#### desde 1986

contagios <- read.csv(text="37,201
                      38,351
                      39,429
                      40,400
                      41,3573
                      42,3958
                      43,2988
                      44,719
                      45,998
                      46,611
                      47,360
                      48,416
                      49,351", header=FALSE)

names(contagios)[1] <- "semana"
names(contagios)[2] <- "infectados"

contagios$semana = contagios$semana

######### El tag de semana es porque el script estaba pensado inicialmente para trabajar con los datos de SIVIGILA
######### Pero para el caso de estos datos, el dato es anual y no semanal... En lugar de semana, debería decir year




############## Dependiendo si los datos son anuales, mensuales, o semanales, se hace una transformacion a días
############## Como en este caso los datos son anuales la transformación es la siguiente:

contagiosxdia <- array(numeric(),c(length(contagios$semana)*365,2)) 
contagiosxdia = data.frame(contagiosxdia)
names(contagiosxdia)[1] <- "dia"
names(contagiosxdia)[2] <- "infectados"

#### Aqui se convierten los datos de años a días ###

for (k in 0:(length(contagios$semana)-1)){
  for (i in 0:364){
    contagiosxdia$infectados[(1+k*365):(1+k*365+i)] = floor((13/365)*(contagios$infectados[k+1]))
    contagiosxdia$dia[(1+k*365+i)] = 365*contagios$semana[k+1]+i-365
  } 
}
contagiosxdia$dia=contagiosxdia$dia+1

#contagiosxdia$dia = contagiosxdia$dia + 13505

#newd=contagiosxdia[order(contagiosxdia$dia, -rank(contagiosxdia$dia), decreasing = FALSE),]

#contagiosxdia = newd

##para datos desde el 2007
################ MODELO SIR ############################
#### Aquí se define el modelo SIR ###

sirvmod = function(t, y, parms) {
  # Pull state variables from y vector
  S1 = y[1]
  S2 = y[2]
  I = y[3]
  R = y[4]
  V = y[5]
  W = y[6]
  C = y[7]
  
  
  ev1 = parms["ev1"] 
  pi1 = parms["pi1"]
  gamma2 = parms["gamma2"] # tasa de recuperacion (sigma en paper)
  N = parms["N"] # total poblacion
  beta2 = parms["beta2"] # es k en el paper
  lambda2 = parms["lambda2"]
  of = parms['of'] # 
  ii = parms['ii'] #
  vr = parms['vr']
  sd = parms['sd']
  
  # social distancing flag
  
  pdis = (W>18330)*(W<(18330+sd*365))#marzo de 2020 al 
  pvac = (W>18330)*(W<(18330+sd*365))#marzo de 2020 al
  
  
  dist <- (W>18330)*(W<18695)
  
  

  brt <- ((W*-1.238913e-06 + 3.566146e-02))/365 # funcion de crecimiento Leandro, birth model
  mu <- ((0.008251563  +  -4.65767e-07*W + -1.9377e-10 *W^2 + 7.939302e-14*W^3 + -1.160313e-17*W^4 + 8.236164e-22*W^5 + -2.837473e-26*W^6 + 3.794165e-31*W^7))/365

  cv1= 0.92*(1-exp(-0.07*(W)/365))*(W>11*365)*(1-vr*pvac)
  
  # modeling variability on the transmision rate
  betaD <- ( beta2*cos(2*pi*W/of) + beta2 )*(1-0.4*pdis)
  
  # total population
  P <- S1 + S2 + I + R + V
  
  
  # Define equations - pi1 PERDIDA DE INMUNIDAD 1
  
  dS1 <- brt*P - betaD*S1*I/P - (1/60)*S1 - mu*S1
  dS2 <- (1-cv1*ev1)*(1/60)*S1 - betaD*S2*I/P + pi1*(V+R) - lambda2*cv1*S2 - mu*S2
  dI <- betaD*(S1+S2)*I/P - gamma2 * I - mu*I
  dR <- gamma2 * I - pi1*R - mu*R
  dV <- (cv1*ev1)*(1/60)*S1 + lambda2*cv1*S2 - pi1*V - mu*V 
  
  dW <-  1
  #dP <- brt*P - mu*P
  dC <- betaD*(S1+S2)*I/P # new cases
  
  res = c(dS1, dS2, dI, dR, dV, dW, dC)
  # Return list of gradients
  list(res)
}


Inis = 4
sirvforc = function(social_dist, vaccination_rate) {
  

  prev1 <- 0.85
  prpi1 <- 1/3650
  prgamma2=0.0011960231
  prbeta2 = 0.0019162239
  prlambda2 = 0.0002000905
  prof=365
  prii=0/365
  
  iniinfect = 50
  pop = 21480065
  iniI = iniinfect
  iniR = 0
  iniE = 0
  iniS = pop-iniI-iniR
  iniW = 1
  iniP = pop
  iniC = 0
  
  
  times = seq(0, 29125, by = 1) 
  parms = c(N = 1, beta2 = prbeta2, lambda2 = prlambda2, gamma2 = prgamma2, ev1=prev1,pi1=prpi1, of=prof, ii=prii, sd = social_dist, vr=vaccination_rate) 
  start = c(S1 = 0, S2 = iniS, I = iniI, R = iniR, V=0, W = iniW, C=iniC)
  

  out=ode(y=start, times=times, func=sirvmod, parms= parms)
  out=as.data.frame(out)
  
  return = out
  
}

##### Acumulacion anual

anualizar = function(vect1) {
  an=replicate(89,0)
  for (i in 1:89){
    an[i] = vect1[365*i]-vect1[365*(i-1)+1]
  }
  return = an
}


outop= sirvforc(1, 0.5)
plot(main="Modelo SIR Pertussis", x=(1:89), y=anualizar(outop$C), ylab="Contagios", xlab= "Year", type="l",col="gray",ylim=c(0,6000))
legend("topright",
       c("datos modelo"),
       fill=c("gray","black")
)

