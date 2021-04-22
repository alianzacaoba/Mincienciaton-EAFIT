#################################
# libraries
#################################

library(deSolve)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())


#############################
# loading measles data
#############################

data.ori <- data.table::fread(input = "population-birth_death_rates-vaccov-1960-2019.csv",
                              header=TRUE, sep=",", na.strings = "",
                              col.names = c("year","population","birth.rate","death.rate",
                                            "vaccov1","vaccov2","cases"))
data <- filter(data.ori, year >= 1980 & year <= 2010)


#################################
# model calibration (1 parameter)
#################################

# measles model
measles_model <- function(time, y, parms){
  
  with(as.list(c(y, parms)),{
    
    # social distancing flag
    #dis.flag <- (W>18330)*(W<18695)
    dis.flag <- 0
    
    # modeling birth rate
    br <- (-1.238913e-06*W + 3.566146e-02)/365
    
    # modeling death rate
    mr <- (8.251563e-03 - 4.657670e-07*W - 1.937700e-10*W^2 + 7.939302e-14*W^3 - 1.160313e-17*W^4 + 8.236164e-22*W^5 - 2.837473e-26*W^6 + 3.794165e-31*W^7)/365
    
    # modeling coverage of the 1st vaccine dose
    v1c <- (W < 3651)*0 + 
      (W >= 3651)*(0.95*(1-exp(-0.1408442*(W-3284)/365)))*(1-0.5*dis.flag)
    #v1c <- 0
    
    # modeling coverage of the 2nd vaccine dose
    v2c <- (W < 9856)*0 + 
      (W >= 9856)*(0.85-exp(-0.4*(W-9489)/365))*(1-0.5*dis.flag)
    #v2c <- 0
    
    # modeling variability on the transmision rate
    beta <- b#( b*cos(2*pi*W/of) + b )*(1-0.4*dis.flag)
    
    # modeling imported cases
    #ii.flag <- (W>7301)*(W<7665)
    #iii <- ii*ii.flag
    
    # adding problematic cases between 1980 to 1995
    #im1 <- ii/365*(W>3650)*(W<8395)
    im1 <- ii/365*(W>3650)*(W<6205)
    im2 <- iii/365*(W>6205)*(W<8760)
    im3 <- iiii/365*(W>17520)
    #im1 <- 0
    #im2 <- 0
    
    # total population
    N <- S + E + I + R + V1 + FF + V2
    
    # states
    dS <- br*N - beta*S*I/N - br*N*v1c - mr*S
    dE <- beta*(S + FF)*I/N - sigma*E - mr*E
    dI <- sigma*E - gamma*I - mr*I + im1 + im2 + im3
    dR <- gamma*I - mr*R
    dV1 <- br*N*v1c*v1e - v2c*V1/(5*365) - mr*V1
    dF <- br*N*v1c*(1 - v1e) - beta*FF*I/N + v2c*(1 - v2e)*V1/(5*365) - mr*FF
    dV2 <- v2c*v2e*V1/(5*365) - mr*V2
    dW <- 1 # counter
    dC <- beta*(S + FF)*I/N + im1 + im2 + im3# new measles cases
    
    # output list
    list(c(dS,dE,dI,dR,dV1,dF,dV2,dW,dC))
  })
}

# evaluation function
sirvforc1 = function(optbeta) {
  
  # parameters
  parameters <- c(b = optbeta, sigma = 1/7, gamma = 1/7,
                  v1e = 0.93, v2e = 0.97,
                  ii=1000, iii=2000, iiii=100)  
  
  # initial conditions
  # Colombia 1980: 9222/26900506*100e3 = 34.3 cases per 100k inhabitants
  init.pob <- 21480065#
  init.infect <- 500
  y <- c(S=init.pob-init.infect, E=0, I=init.infect, R=0, V1=0, FF=0, V2=0, W=0, C=0)
  
  # time frame
  start.year <- 1970
  final.year <- 2010
  number.years <- final.year - start.year + 1
  end.time <- number.years*365
  time <- 1:end.time
  
  # solving model
  model.sol <- ode(y, time, measles_model, parameters, method = "ode2")
  
  # first and last day of each year from 1970 to 2029
  start.end.annual <- lapply(0:(number.years-1), function(x) c(1,365)+365*x)
  start.end.annual <- unlist(start.end.annual)
  
  # cases per year
  model.sol.2 <- model.sol %>% data.frame() %>%
    mutate(year = rep(start.year:final.year, each=365)) %>%
    filter(time %in% start.end.annual) %>%
    group_by(year) %>% 
    summarise(total.cases = abs(diff(C)))
  
  # returning total cases between 1994 & 2010
  return(model.sol.2 %>% filter(year >= 1980) %>% pull(total.cases))
}

# calibration function
LSmodelo1 = function(optbeta) {
  loss = sqrt(sum((sirvforc1(optbeta)-data$cases)^2))
  return(loss)
}

# performing calibration
optbeta <- 0.15
calib.out.1 <- optim(par=optbeta, LSmodelo1, method="Brent", lower = 0.10, upper=0.2)


#################################
# visualizing results
#################################

# parameters
params <- c(b = calib.out.1$par, # last optim: 0.162596
            sigma = 1/7, gamma = 1/7,
            #v1c = 0.95, v2c = 0.95,
            v1e = 0.93, v2e = 0.97,
            ii=1000, iii=2000, iiii=100)

# initial conditions
# Colombia 1980: 9222/26900506*100e3 = 34.3 cases per 100k inhabitants
init.pob <- 21480065#
init.infect <- 500
y <- c(S=init.pob-init.infect, E=0, I=init.infect, R=0, V1=0, FF=0, V2=0, W=0, C=0)

# time frame
start.year <- 1970
final.year <- 2020
number.years <- final.year - start.year + 1
end.time <- number.years*365
time <- 1:end.time

# solving model
model.sol <- ode(y, time, measles_model, params, method = "ode2") 

# 
# first and last day of each year from 1970 to 2029
start.end.annual <- lapply(0:(number.years-1), function(x) c(1,365)+365*x)
start.end.annual <- unlist(start.end.annual)
model.sol.2 <- model.sol %>% data.frame() %>%
  mutate(year = rep(start.year:final.year, each=365)) %>%
  filter(time %in% start.end.annual) %>%
  group_by(year) %>% 
  summarise(cases = abs(diff(C)))

# building dataframe for ggplot
df.plot <- rbind.data.frame(filter(data.ori, year >= 1980 & year <= 2019) %>%
                              select(year,cases) %>% mutate(trend = "True"),
                            model.sol.2 %>% filter(year >= 1980 & year <= 2019) %>%
                              mutate(trend = "Simulated"))

#
base.plot <- df.plot %>% ggplot() +
  geom_line(aes(x=year, y=cases, color=trend), size = 1) +
  geom_point(aes(x=year, y=cases, color=trend), size = 2) +
  theme(text = element_text(size=12),
        legend.title = element_blank()) +
  xlab("Years") + ylab("Total cases") +
  scale_colour_manual(labels = c("Simulated cases", "Reported cases"),
                      values = c("black",rgb(red=1, green=0, blue=0, alpha=0.5)))
  #theme(legend.position="bottom") +
  

####################
# report figure
####################

# calling pdf printer
pdf(file = "meales-1-calibration.pdf", width = 6, height = 4)

base.plot +
  theme(legend.position = c(0.855, 0.895),
        legend.background = element_rect(size=0.3, linetype="solid", colour ="black"))

# turn off the pdf printer
dev.off()
