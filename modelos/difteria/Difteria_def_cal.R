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
  
  # social distancing flag
  dist <- (W>18330)*(W<18695)
  #dist <- 0
  
  # Fixed variables modeling (Based on population-birth_death_rates-vaccov-1960-2019 and meales_data_fit_models)
  # Birth rate
  brt <- ((W*-1.238913e-06 + 3.566146e-02))/365 # funcion de crecimiento Leandro, birth model
  # Death rate
  mu <- ((0.008251563  +  -4.65767e-07*W + -1.9377e-10 *W^2 + 7.939302e-14*W^3 + -1.160313e-17*W^4 + 8.236164e-22*W^5 + -2.837473e-26*W^6 + 3.794165e-31*W^7))/365
  # dose 1 vaccine coverage (base coverage input 1)
  cv1 <- 0.92*(1-exp(-0.07*(W)/365))*(W>11*365)#*(1-vr*pvac)

  
  # modeling variability on the transmision rate
  betaD <- ( beta2*cos(2*pi*W/of) + beta2 )*(1-0.4*dist)
  
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
#                    calibration function                                      #
#                       for Beta and ii                                        #
################################################################################
calib_bii_function <- function(optbeta, optii) {
  
  
  # Parameters:
  prev1 <- 1 # Efectivity of dose 1
  prgamma2=1/14 # Tasa de recuperación 
  prbeta2 = abs(optbeta) # beta
  
  # Initial state: y
  iniinfect = 226 # infectados iniciales (1997)
  pop = 21480065 # población de 1997
  iniI = iniinfect 
  iniR = 0
  iniS = iniI-iniR
  iniE = 0
  iniW = 0
  ii = abs(optii)
  
  parms = c(beta2 = prbeta2, gamma2 = prgamma2, 
            ev1=prev1, 
            of = 365, 
            ii=ii,eta_v=1/(365*30),eta_R=1/(365*30)) 
  start = c(S = pop, I = iniI, R = iniR, V=0, W = iniW, C=0)
  
  
  # time frame
  start.year <- 1970
  final.year <- 2018
  number.years <- final.year - start.year + 1
  end.time <- number.years*365
  time <- 1:end.time
  
  #################################
  # solving model
  #################################
  
  out <- ode(y=start, times=time, func=Difteria_30, 
             parms= parms)
  out <- as.data.frame(out) # model output as dataframe.
  
  # first and last day of each year from 1970 to 2029
  start.end.annual <- lapply(0:(number.years-1), function(x) c(1,365)+365*x)
  start.end.annual <- unlist(start.end.annual)
  
  # cases per year
  out1 <- out %>% data.frame() %>%
    mutate(year = rep(start.year:final.year, each=365)) %>%
    filter(time %in% start.end.annual) %>%
    group_by(year) %>% 
    summarise(total.cases = abs(diff(C)))
  
  # returning total cases between 1994 & 2010
  return(out1 %>% filter(year >= 1974) %>% pull(total.cases))
  
  
}

################################################################################
#                    calibration function                                      #
#                       for Beta                                               #
################################################################################
calib_b_function <- function(optbeta) {
  
  
  # Parameters:
  prev1 <- 1 # Efectivity of dose 1
  prgamma2=1/14 # Tasa de recuperación 
  prbeta2 = abs(optbeta) # beta
  
  # Initial state: y
  iniinfect = 226 # infectados iniciales (1997)
  pop = 21480065 # población de 1997
  iniI = iniinfect 
  iniR = 0
  iniS = iniI-iniR
  iniE = 0
  iniW = 0
  ii = 0.003
  
  parms = c(beta2 = prbeta2, gamma2 = prgamma2, 
            ev1=prev1, 
            of = 365, 
            ii=ii,eta_v=1/(365*30),eta_R=1/(365*30)) 
  start = c(S = pop, I = iniI, R = iniR, V=0, W = iniW, C=0)
  
  
  # time frame
  start.year <- 1970
  final.year <- 2018
  number.years <- final.year - start.year + 1
  end.time <- number.years*365
  time <- 1:end.time
  
  #################################
  # solving model
  #################################
  
  out <- ode(y=start, times=time, func=Difteria_30, 
             parms= parms)
  out <- as.data.frame(out) # model output as dataframe.
  
  # first and last day of each year from 1970 to 2029
  start.end.annual <- lapply(0:(number.years-1), function(x) c(1,365)+365*x)
  start.end.annual <- unlist(start.end.annual)
  
  # cases per year
  out1 <- out %>% data.frame() %>%
    mutate(year = rep(start.year:final.year, each=365)) %>%
    filter(time %in% start.end.annual) %>%
    group_by(year) %>% 
    summarise(total.cases = abs(diff(C)))
  
  # returning total cases between 1994 & 2010
  return(out1 %>% filter(year >= 1974) %>% pull(total.cases))
  
  
}


################################################################################
#                     To call calibration function                             #
# There are two possibilities, calibrate both beta and ii or just calibrate    #
# beta. There are two loss functions, Least Squares (LS) or Maximum Likelihood #
# Estimator (MLE).                                                             #
# First thing to do is load real cases data.                                   #
# Optimal parameters will be stored on calib_beta and calib_ii                 #
################################################################################

############################
# Loading real cases data  #
# Data taken from OMS      #
############################

data <- read.csv("C:/Users/sulon/OneDrive - Universidad EAFIT/modelo-difteria/cases_per_year_diphtheria_OMS.csv", 
                 sep=";", fileEncoding="UTF-8-BOM")

# Data from 1990
#data <- read.csv("C:/Users/sulon/OneDrive - Universidad EAFIT/modelo-difteria/diphtheria_since_90.csv", 
#                                  sep=";", fileEncoding="UTF-8-BOM")

################################################################################
# For beta and imported infectious. (beta and ii)

##################################
# With LS model,                 #
#       for beta and ii          #
##################################

LSmodelo1 = function(vect1) {
  loss = sqrt(sum((calib_bii_function(vect1[1],vect1[2])-data$cases)^2))
  return(loss)
}

# To run this, uncomment the following lines: 
#optbeta <- 0.01
#optii <- 0.001
#calib_out1 <- optim(par=c(1,1),fn=LSmodelo1,method = "L-BFGS-B",lower = c(0,0), upper = c(1,0.002))
#calib_beta <- calib_out1$par[1]
#calib_ii <- calib_out1$par[2]


#################################
#  With MLENorm,                #
#         for beta and ii       #
#################################

MLENormmodelo = function(vect1) { # Maximun likelihood estimation
   loss = -(sum(dnorm(data$cases,mean=calib_bii_function(vect1[1],vect1[2]),sd=100,log=TRUE)))
   return = loss
}

# To run this, uncomment the following lines: 
optbeta <- 0.0594873493117294
optii <- 0.003 
calib_out1 <- optim(par=c(1,1),fn=MLENormmodelo,method = "L-BFGS-B",lower = c(0,0), upper = c(1,0.003))
calib_beta <- calib_out1$par[1]
calib_ii <- calib_out1$par[2]

#
################################################################################
################################################################################
# For beta.
#################################
# With LS model                 #
#         for beta              #
#################################

LSmodelo1 = function(optbeta) {
  loss = sqrt(sum((calib_b_function_function(vect1[1])-data$cases)^2))
  return(loss)
}

# To run this, uncomment the following lines: 
#optbeta <- 0.01
#calib_out1 <- optim(par= optbeta,fn=LSmodelo1,method = "L-BFGS-B",lower = 0, upper = 2)
#calib_beta <- calib_out1$par

#################################
# With MLENorm                  #
#         for beta              #
#################################

MLENormmodelo = function(optbete) { # Maximun likelihood estimation
  loss = -(sum(dnorm(data$cases,mean=calib_b_function(optbete),sd=100,log=TRUE)))
  return = loss
}

# To run this, uncomment the following lines: 
#optbeta <- 0.01
#calib_out1 <- optim(par=optbeta , fn=MLENormmodelo, method = "L-BFGS-B",lower = 0, upper = 1)
#calib_beta <- calib_out1$par

#
################################################################################

################################################################################
# Running the model with optimal beta and ii parameters                        #
################################################################################
####### Con MLE
#### Beta calibrado desde 1974 -> 0.0594873493117294
#### ii calibrado desde 1974 -> 0.003

####### Con MLE
#### Beta calibrado desde 1990 -> 0.0315473125807805
#### ii calibrado desde 1990 -> 0.003

####### Con MLE
#### Beta calibrado desde 1990 -> 0.0485222825597622
#### ii calibrado desde 1990 -> 0.003
calib_beta <- 0.0594873493117294
calib_ii <- 0.003
# Parameters:
prev1 <- 1 # Efectivity of dose 1
prgamma2=1/14 # Tasa de recuperación 
prbeta2 = calib_beta # beta 
# Initial state: y
iniinfect = 226 # infectados iniciales (1997)
pop = 21480065 # población de 1997
iniI = iniinfect 
iniR = 0
iniS = iniI-iniR
iniW = 0
ii = calib_ii

parms = c(beta2 = prbeta2, gamma2 = prgamma2, 
          ev1=prev1, 
           of = 365, 
          ii=ii,eta_v=1/(365*30),eta_R=1/(365*30)) 
start = c(S = pop, I = iniI, R = iniR, V=0, W = iniW, C=0)


# time frame
start.year <- 1970
final.year <- 2019
number.years <- final.year - start.year + 1
end.time <- number.years*365
time <- 1:end.time

#################################
# solving model
#################################

model.sol = ode(y=start, times=time, func=Difteria_30, parms= parms)
out = as.data.frame(model.sol)


casos_por_100k <-NULL # c?lcula base por a?os de casos por 100.000 habitantes y casos totales
for (i in 1:number.years) { 
  fin = i*365
  #Casos
  tpop = out[fin,2]+out[fin,3]+out[fin,4]+out[fin,5]
  cas = (out$C[i*365]-out$C[i*365-364])*100000/tpop
  cas0 = (out$C[i*365]-out$C[i*365-364])
  ## Totales
  rec = data.frame(ano=i+1969,casos_ci=cas,casos_t=cas0)
  casos_por_100k<-rbind(casos_por_100k,rec)
}



#################################
# visualization of results
#################################

library(ggplot2)
casos_observados <- casos_por_100k %>% filter(ano >= 1974 & ano <= 2018)
casos_reales <- data$cases
color_group <- c("black","red")
x_ax = (1974:2018)

ggplot() +
  geom_line(aes(x=x_ax,y=casos_observados$casos_t, color='Predicted cases')) + geom_point(aes(x=x_ax,y=casos_observados$casos_t, color='Predicted cases')) +
  geom_line(aes(x=x_ax,y=casos_reales, color='Observed cases')) + geom_point(aes(x=x_ax,y=casos_reales, color='Observed cases')) +
  scale_colour_manual(values=color_group) +
  labs(y='Anual Cases', x = "Years",colour="") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.8)) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(hjust=1)) + 
  ggtitle("Diptheria Model Fitting")



datos<-(c(rep(NA,4),casos_reales, rep(NA,1)))

BASE_OUT<-data.frame(Unidades_temporales=casos_por_100k$ano,
                     id_variable=rep("UE-DIF-C",50),
                     Valor_proyectado= casos_por_100k$casos_t,
                     "Proyectado/100.000 hab"=casos_por_100k$casos_ci,
                     "Valor observado"=datos,
                     Error=abs(casos_por_100k$casos_t-datos)/datos)

library("xlsx")
write.xlsx(BASE_OUT, "UE_DIF_OBS_PRO.xlsx")





# 
# # plotting only infected people
# plot(out[,"I"], col = "red", type = "l")
