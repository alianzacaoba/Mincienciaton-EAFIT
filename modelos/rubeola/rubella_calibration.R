# calibracion MLE 

library(deSolve)
library(dfoptim)
options(scipen=999)

############### Aquí se cargan los datos de contagios ############################
#### Los datos que se suministran aquí son datos de prueba para la calibracion, 
#### tomados del numero de reportes de casos por a�os de la OMS.

contagios <- read.csv(text="27,226
                      28,1906
                      29,974
                      30,157
                      31,70
                      32,1206
                      33,47
                      34,45
                      35,85
                      36,6
                      37,2
                      38,4
                      39,4
                      40,0
                      41,1
                      42,1
                      43,0
                      44,0
                      45,0
                      46,0
                      47,0
                      48,0
                      49,0", header=FALSE)

names(contagios)[1] <- "semana"
names(contagios)[2] <- "infectados"

contagios$semana = contagios$semana

# El tag de semana es porque el script estaba pensado inicialmente para
# trabajar con los datos de SIVIGILA Pero para el caso de estos datos, 
# el dato es anual y no semanal... En lugar de semana, debería decir year


# Dependiendo si los datos son anuales, mensuales, o semanales, se hace una 
# transformacion a días Como en este caso los datos son anuales 
# la transformación es la siguiente:

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
  
  # social distancing flag
  #dist <- (W>18330)*(W<18695)
  dist <- 0
  
  # Fixed variables modeling (Based on population-birth_death_rates-vaccov-1960-2019 and meales_data_fit_models)
  # Birth rate
  brt <- ((W*-1.238913e-06 + 3.566146e-02))/365 # funcion de crecimiento Leandro, birth model
  # Death rate
  mu <- ((0.008251563  +  -4.65767e-07*W + -1.9377e-10 *W^2 + 7.939302e-14*W^3 + -1.160313e-17*W^4 + 8.236164e-22*W^5 + -2.837473e-26*W^6 + 3.794165e-31*W^7))/365
  # dose 1 vaccine coverage (base coverage input 1)
  cv1 <- (W < 3651)*0 + (W >= 3651)*(0.95*(1-exp(-0.1408442*(W-3284)/365)))*(1-0.5*dist)
  #cv1=0

  # dose 2 vaccine coverage
  cv2 <- (W < 9856)*0 + (W >= 9856)*(0.85-exp(-0.4*(W-9489)/365))*(1-0.5*dist)
  #cv2=0
  
  # modeling variability on the transmision rate
  betaD <- ( beta2*cos(2*pi*W/of) + beta2 )*(1-0.4*dist)
  
  # total population
  P <- S + E + I + R + V + V2
  
  # Define equations
  dS <- brt*(1-cv1*ev1)*P + V*(1/(365*4))*cv2*(1-ev2) - betaD*S*I/P - ee*brt*E - ie*brt*E - mu*S #br*N - beta*S*I/N - br*N*v1c - mr*S
  dE <- betaD*S*I/P + ee*brt*E + ie*brt*E - e*E - mu*E
  dI <- e*E - (mu + gamma2)*I + ii
  dR <- gamma2 * I - mu*R
  dV <- brt*cv1*ev1*P - V*(1/(365*4))*cv2 - mu*V
  dV2 <- V*(1/(365*4))*cv2*ev2 - mu*V2
  dW <- 1 
  #dP <- brt*P - mu*P
  dC <- e*E + ii # new measles cases
  
  res = c(dS, dE, dI, dR, dV, dV2, dW, dC)
  # Return list of gradients
  list(res)
}


############################################################################################################
#### Aquí se pone la función que corre el modelo, de forma que los argumentos sean beta y gamma y el resto de los
#### parametros sean estáticos
#### notese en la linea times = seq(0, contagiosxdia$dia[length(contagiosxdia$dia)], by = 1)
#### que el modelo está hecho para que corra la misma cantidad de dias disponibles en los datos, y el by = 1
#### es importante ya que el script esta hecho para tomar directamente el out de infectados como un dia
#### en cada valor del data frame
Inis = 4
sirvforc = function(optbeta,optgamma) {
  
  prmu = 1/((50000000/277177)*365) 
  #prv=1/((50000000/637669)*365)
  #coberturav <- input$num1 * (1/100) # 0 - 0, 1 - 0 ... 3651 - 5%; 3652 - 5.004%, 17885 - 90% 
  prev1 <- 0.88
  prev2 <- 0.97
  #prp = eficaciav*coberturav
  prgamma2=abs(optgamma)
  # N 
  prbeta2 = abs(optbeta)
  pree = 0.65 # de acuerdo al paper
  prie = 0.65 # de acuerdo al paper
  pre = 1/90 # de acuerdo al paper
  prof=365
  prii=0/365
  
  iniinfect = 5
  pop = 21480065
  iniI = iniinfect
  iniR = 0
  iniE = 0
  iniS = pop-iniI-iniR
  iniW = 1
  iniP = pop
  iniC = 0
  #pop <- W*(1642.828) + 21184610.465
  
  times = seq(0, 19125, by = 1) # aca esta hasta 31 de dic de 2021
  parms = c(mu = prmu, N = 1, beta2 = prbeta2, gamma2 = prgamma2, ev1=prev1, ev2 = prev2, ee = 0.65, ie = 0.65, e = pre, of=prof, ii=prii) 
  start = c(S = iniS, E = iniE, I = iniI, R = iniR, V=0, V2=0, W = iniW, C=iniC)
  
  #parms2 = c(mu = prmu, N = 1, beta2 = prbeta2, gamma2 = prgamma2, ev1=prev1, ev2 = prev2, ee = 0.65, ie = 0.65, e = pre) 
  #start2 = c(S = iniS, E = iniE, I = iniI, R = iniR, V=0, V2=0)
  
  out=ode(y=start, times=times, func=sirvmod, parms= parms)
  out=as.data.frame(out)
  
  return = out
  
}

##### Acumulacion anual

anualizar = function(vect1) {
  an=replicate(49,0)
  for (i in 1:49){
    an[i] = vect1[365*i]-vect1[365*(i-1)+1]
  }
  return = an
}


##### Aqui se hacen las funciones para poner en optim()

LSmodelo = function(vect1) {
  loss = sqrt((sum(((anualizar(sirvforc(vect1[1],vect1[2])$I)[27:49])-((contagios$infectados)))))^2)
  return = loss
}

MLENormmodelo = function(vect1) {
  loss = -(sum(dnorm(contagiosxdia$infectados,mean=sirvforc(vect1[1],vect1[2])$I[9491:(9491+(length(contagiosxdia$infectados)-1))],sd=100,log=TRUE)))
  return = loss
}

MLEPoissmodelo = function(vect1) {
  loss = -(sum(dpois(contagiosxdia$infectados,lambda = sirvforc(vect1[1],vect1[2])$I[9491:(9491+(length(contagiosxdia$infectados)-1))] ,log=TRUE)))
  return = loss
}
MLEPoissmodelo3 = function(vect1) {
  loss = -(sum(dpois(contagios$infectados,lambda = (anualizar(sirvforc(vect1[1],vect1[2])$C)[(27:49)]) ,log=TRUE)))
  return = loss
}

### A continuacion tenemos las calibraciones para dos paràmetros
optbeta=100
opte=100
optgamma=100
#optim(par=c(1,1),fn=LSmodelo, method="Nelder-Mead")
#optim(par=c(1,1), MLENormmodelo,  method = "Nelder-Mead")
#optim(par=c(1,1), fn=MLEPoissmodelo,  method="Nelder-Mead")
calib_out = optim(par=c(1,1), fn=MLEPoissmodelo3,  method="Nelder-Mead")

### Aquí graficamos manualmente con uno de los outputs, esta parte es solo para hacer seguimiento de los resultados:


#### Grafico con LS

#outop = sirvforc(0.0332,0.0318,5322.05)
#outop= sirvforc(82.18, 435.262)
#plot(main="Modelo SIR con demografía y vacunación VS datos reales", x=outop$time, y=outop$I, ylab="Contagios", xlab= "Días", type="l",col="gray", xlim = c(0,12000), ylim=c(0,80))
#lines(x=contagiosxdia$dia, y=contagiosxdia$infectados, type = "l")
#legend("topright",
#       c("datos modelo","datos reales"),
#       fill=c("gray","black")
#)

#### Grafico con LS

#outop = sirvforc(0.0332,0.0318,5322.05)
#outop= sirvforc(4263, 3595)
#plot(main="Modelo SIR con demografía y vacunación VS datos reales", x=outop$time, y=outop$I, ylab="Contagios", xlab= "Días", type="l",col="gray", xlim = c(0,12000), ylim=c(0,80))
#lines(x=contagiosxdia$dia, y=contagiosxdia$infectados, type = "l")
#legend("topright",
#       c("datos modelo","datos reales"),
#       fill=c("gray","black")
#)

#outop = sirvforc(0.0332,0.0318,5322.05)
# beta_optimo = 0.01681578 
# gamma_optimo = 0.01340360
outop= sirvforc(0.01681578, 0.01340360)
plot(main="Modelo SIR con demografía y vacunación VS datos reales", x=(1:49), y=anualizar(outop$C), ylab="Contagios", xlab= "Year", type="l",col="gray",ylim=c(0,4000))
points(x=contagios$semana, y=contagios$infectados, type = "p")
legend("topright",
       c("datos modelo","datos reales"),
       fill=c("gray","black")
)


plot(main="Modelo SIR con demografía y vacunación VS datos reales", x=outop$time, y=outop$I, ylab="Contagios", xlab= "Días", type="l",col="gray")

(c(outop$W[17521], outop$S[17521], outop$V[17521], outop$V2[17521], outop$E[17521], outop$I[17521], outop$R[17521]))


library(ggplot2)
casos_observados <- data.frame(casos_obs = anualizar(outop$C)[26:48])
casos_reales <- data.frame(casos_re = contagios$infectados)
color_group <- c("black","red")

ggplot() +
  geom_line(aes(x=(1997:2019),y=casos_observados$casos_obs, color='Predicted cases')) + geom_point(aes(x=(1997:2019),y=casos_observados$casos_obs, color='Predicted cases')) +
  geom_line(aes(x=(1997:2019),y=casos_reales$casos_re, color='Observed cases')) + geom_point(aes(x=(1997:2019),y=casos_reales$casos_re, color='Observed cases')) +
  scale_colour_manual(values=color_group) +
  labs(y='Anual Cases', x = "Years",colour="") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.8)) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1)) + 
  ggtitle("Rubella Model Fitting")


################################################################################
# Salida de los datos                                                          #
################################################################################
casos_observados <- data.frame(casos_obs = anualizar(outop$C))
casos_reales <- data.frame(casos_re = contagios$infectados)

casos_por_100k <-NULL # c?lcula base por a?os de casos por 100.000 habitantes y casos totales
for (i in 1:50) { 
  fin = i*365
  #Casos
  tpop = outop[fin,2]+outop[fin,3]+outop[fin,4]+outop[fin,5]+outop[fin,6]+outop[fin,7]
  cas = (outop$C[i*365]-outop$C[i*365-364])*100000/tpop
  cas0 = (outop$C[i*365]-outop$C[i*365-364])
  ## Totales
  rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0)
  casos_por_100k<-rbind(casos_por_100k,rec)
}


datos<-(c(rep(NA,26),casos_reales$casos_re, rep(NA,1)))

BASE_OUT<-data.frame(Unidades_temporales=casos_por_100k$ano,
                     id_variable=rep("UE-RUB-C",50),
                     Valor_proyectado= casos_por_100k$casos_t,
                     "Proyectado/100.000 hab"=casos_por_100k$casos_ci,
                     "Valor observado"=datos,
                     Error=abs(casos_por_100k$casos_t-datos)/datos)

library("xlsx")
write.xlsx(BASE_OUT, "Rub_Observado-Pronostico.xlsx")



