#################################
# libraries
#################################

library(deSolve)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())


#################################
# base measles model
#################################

measles_model_base <- function(time, y, parms){
  
  with(as.list(c(y, parms)),{
    
    # social distancing flag
    # (also used for setting periods with decreased vaccine coverage)
    dis.flag <- (W>18330)*(W<(18330+tdis*365))

    # modeling birth rate
    br <- (-1.238913e-06*W + 3.566146e-02)/365
    
    # modeling death rate
    mr <- (8.251563e-03 - 4.657670e-07*W - 1.937700e-10*W^2 + 7.939302e-14*W^3 - 1.160313e-17*W^4 + 8.236164e-22*W^5 - 2.837473e-26*W^6 + 3.794165e-31*W^7)/365
    
    # modeling coverage of the 1st vaccine dose
    v1c <- (W < 3651)*0 + 
      (W >= 3651)*(0.95*(1-exp(-0.1408442*(W-3284)/365)))*(1-rvac*dis.flag)
    #v1c <- 0
    
    # modeling coverage of the 2nd vaccine dose
    v2c <- (W < 9856)*0 + 
      (W >= 9856)*(0.85-exp(-0.4*(W-9489)/365))*(1-rvac*dis.flag)
    #v2c <- 0
    
    # modeling variability on the transmision rate
    beta <- b*(1-0.4*dis.flag)#( b*cos(2*pi*W/of) + b )*(1-0.4*dis.flag)
    

    # adding problematic cases
    im1 <- ii/365*(W>3650)*(W<6205)
    im2 <- iii/365*(W>6205)*(W<8760)
    im3 <- iiii/365*(W>17520)#*(W<19345)
    
    
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
    dC <- beta*(S + FF)*I/N + im1 + im2 + im3 # new measles cases
    
    # output list
    list(c(dS,dE,dI,dR,dV1,dF,dV2,dW,dC))
  })
}

#################################
# encapsulated model
#################################

# tdis = social distance time frame (portion of a year, between 0 and 1)
# rvac = reduction of vaccine coverage (between 0 and 1)
measles_model_cases <- function(tdis,rvac){
  
  # parameters
  parameters <- c(b = 0.162596, sigma = 1/7, gamma = 1/7,
                  #v1c = 0.95, v2c = 0.95,
                  v1e = 0.93, v2e = 0.97,
                  ii=1000, iii=2000, iiii=100,
                  tdis=tdis, rvac=rvac) 
  
  # initial conditions
  # Colombia 1980: 9222/26900506*100e3 = 34.3 cases per 100k inhabitants
  init.pob <- 21480065     #26900506#
  init.infect <- 500#
  y <- c(S=init.pob-init.infect, E=0, I=init.infect, R=0, V1=0, FF=0, V2=0, W=0, C=0)
  
  # time frame
  start.year <- 1970
  final.year <- 2029
  number.years <- final.year - start.year + 1
  end.time <- number.years*365
  time <- 1:end.time
  
  # solving model
  model.sol <- ode(y, time, measles_model_base, parameters, method = "ode2")
  
  # total population by year
  total.pop.annual <- model.sol %>% data.frame() %>%
    filter(time %in% seq(365, 365*number.years, 365)) %>%
    select(-c(time,W,C)) %>%
    apply(.,1,sum)
  
  # first and last day of each year from 1970 to 2029
  start.end.annual <- lapply(0:(number.years-1), function(x) c(1,365)+365*x)
  start.end.annual <- unlist(start.end.annual)
  
  # new cases per year (total & per 100k inhabitants)
  model.sol.2 <- model.sol %>% data.frame() %>%
    mutate(year = rep(start.year:final.year, each=365)) %>%
    filter(time %in% start.end.annual) %>%
    group_by(year) %>% 
    summarise(total.cases = abs(diff(C)))
  model.sol.2$total.cases.per.100k <- 100e3*model.sol.2$total.cases/total.pop.annual
  
  # output
  return(
    list(model.sol = model.sol,
         model.sol.2 = model.sol.2)
  )
}


#################################
# running diverse scenarios
#################################

# without covid-19 (baseline)
model.sol.0m.0r <- measles_model_cases(tdis=0, rvac=0)

# half a year with social distancing
# and diverse reduction of vaccine coverage
model.sol.6m.10r <- measles_model_cases(tdis=0.5, rvac=0.1)
model.sol.6m.25r <- measles_model_cases(tdis=0.5, rvac=0.25)
model.sol.6m.50r <- measles_model_cases(tdis=0.5, rvac=0.5)

# a year with social distancing
# and diverse reduction of vaccine coverage
model.sol.12m.10r <- measles_model_cases(tdis=1, rvac=0.1)
model.sol.12m.25r <- measles_model_cases(tdis=1, rvac=0.25)
model.sol.12m.50r <- measles_model_cases(tdis=1, rvac=0.5)



#################################
# organizing results
#################################

# all of these are dataframes that will be used
# by ggplot

df0 <- model.sol.0m.0r$model.sol.2
df1 <- model.sol.6m.10r$model.sol.2
df2 <- model.sol.6m.25r$model.sol.2
df3 <- model.sol.6m.50r$model.sol.2
df4 <- model.sol.12m.10r$model.sol.2
df5 <- model.sol.12m.25r$model.sol.2
df6 <- model.sol.12m.50r$model.sol.2

type.simulations <- c("Six months (coverage -10%)","Six months (coverage -25%)",
                      "Six months (coverage -50%)",
                      "One year (coverage -10%)","One year (coverage -25%)",
                      "One year (coverage -50%)","Baseline")

df0$Simulation <- type.simulations[7]
df1$Simulation <- type.simulations[1]
df2$Simulation <- type.simulations[2]
df3$Simulation <- type.simulations[3]
df4$Simulation <- type.simulations[4]
df5$Simulation <- type.simulations[5]
df6$Simulation <- type.simulations[6]

# final dataframe for ggplot
df <- rbind(df0,df1,df2,df3,df4,df5,df6)
df$Simulation <- factor(df$Simulation, levels = type.simulations)


#####################################
# visualizing results I
#####################################

line.cols <- c("Baseline" = "black",
               "Six months (coverage -10%)" = "blue",
               "Six months (coverage -25%)" = "green",
               "Six months (coverage -50%)" = "red",
               "One year (coverage -10%)" = "blue",
               "One year (coverage -25%)" = "green",
               "One year (coverage -50%)" = "red")
lines.typ <- c("Baseline" = "solid",
               "Six months (coverage -10%)" = "solid",
               "Six months (coverage -25%)" = "solid",
               "Six months (coverage -50%)" = "solid",
               "One year (coverage -10%)" = "dotted",
               "One year (coverage -25%)" = "dotted",
               "One year (coverage -50%)" = "dotted")


# calling pdf printer
#pdf(file = "measles_model_scenarios-3-2018-2029.pdf", width = 6, height = 4)
pdf(file = "measles_model_scenarios-3-2024-2029.pdf", width = 6, height = 4)


# visualazing results
df %>%
  filter(year >= 2024) %>%
  #filter(year >= 2018) %>% #2018, 2024
  ggplot() +
  geom_line(aes(x=year, y=total.cases, color=Simulation,
                linetype=Simulation), size = 1) +
  #ggtitle("Measles cases (2018-2029)") +
  ggtitle("Measles cases (2024-2029)") +
  ylab("Annual cases")  + xlab("") +
  #scale_x_continuous(breaks = seq(2018,2029, by = 1)) +
  scale_x_continuous(breaks = seq(2024,2029, by = 1)) +
  theme(text = element_text(size=10)) +
  scale_color_manual(values = line.cols) +
  scale_linetype_manual(values = lines.typ) + 
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=3))


# turn off the pdf printer
dev.off()


#####################################
# organizing results II
#####################################

daily.total.model.pop <- model.sol.0m.0r$model.sol %>% data.frame() %>%
  select(-c(time,W,C)) %>% rowSums()

df.proport.0.V1 <- data.frame(Proportion=model.sol.0m.0r$model.sol[,"V1"]/daily.total.model.pop,
                              Simulation=type.simulations[7],
                              Vaccine="MMR1",
                              time=model.sol.0m.0r$model.sol[,"time"])
df.proport.0.V2 <- data.frame(Proportion=model.sol.0m.0r$model.sol[,"V2"]/daily.total.model.pop,
                              Simulation=type.simulations[7],
                              Vaccine="MMR2",
                              time=model.sol.0m.0r$model.sol[,"time"])

df.proport.1.V1 <- data.frame(Proportion=model.sol.6m.10r$model.sol[,"V1"]/daily.total.model.pop,
                              Simulation=type.simulations[1],
                              Vaccine="MMR1",
                              time=model.sol.0m.0r$model.sol[,"time"])
df.proport.1.V2 <- data.frame(Proportion=model.sol.6m.10r$model.sol[,"V2"]/daily.total.model.pop,
                              Simulation=type.simulations[1],
                              Vaccine="MMR2",
                              time=model.sol.0m.0r$model.sol[,"time"])

df.proport.2.V1 <- data.frame(Proportion=model.sol.6m.25r$model.sol[,"V1"]/daily.total.model.pop,
                              Simulation=type.simulations[2],
                              Vaccine="MMR1",
                              time=model.sol.0m.0r$model.sol[,"time"])
df.proport.2.V2 <- data.frame(Proportion=model.sol.6m.25r$model.sol[,"V2"]/daily.total.model.pop,
                              Simulation=type.simulations[2],
                              Vaccine="MMR2",
                              time=model.sol.0m.0r$model.sol[,"time"])

df.proport.3.V1 <- data.frame(Proportion=model.sol.6m.50r$model.sol[,"V1"]/daily.total.model.pop,
                              Simulation=type.simulations[3],
                              Vaccine="MMR1",
                              time=model.sol.0m.0r$model.sol[,"time"])
df.proport.3.V2 <- data.frame(Proportion=model.sol.6m.50r$model.sol[,"V2"]/daily.total.model.pop,
                              Simulation=type.simulations[3],
                              Vaccine="MMR2",
                              time=model.sol.0m.0r$model.sol[,"time"])

df.proport.4.V1 <- data.frame(Proportion=model.sol.12m.10r$model.sol[,"V1"]/daily.total.model.pop,
                              Simulation=type.simulations[4],
                              Vaccine="MMR1",
                              time=model.sol.0m.0r$model.sol[,"time"])
df.proport.4.V2 <- data.frame(Proportion=model.sol.12m.10r$model.sol[,"V2"]/daily.total.model.pop,
                              Simulation=type.simulations[4],
                              Vaccine="MMR2",
                              time=model.sol.0m.0r$model.sol[,"time"])

df.proport.5.V1 <- data.frame(Proportion=model.sol.12m.25r$model.sol[,"V1"]/daily.total.model.pop,
                              Simulation=type.simulations[5],
                              Vaccine="MMR1",
                              time=model.sol.0m.0r$model.sol[,"time"])
df.proport.5.V2 <- data.frame(Proportion=model.sol.12m.25r$model.sol[,"V2"]/daily.total.model.pop,
                              Simulation=type.simulations[5],
                              Vaccine="MMR2",
                              time=model.sol.0m.0r$model.sol[,"time"])

df.proport.6.V1 <- data.frame(Proportion=model.sol.12m.50r$model.sol[,"V1"]/daily.total.model.pop,
                              Simulation=type.simulations[6],
                              Vaccine="MMR1",
                              time=model.sol.0m.0r$model.sol[,"time"])
df.proport.6.V2 <- data.frame(Proportion=model.sol.12m.50r$model.sol[,"V2"]/daily.total.model.pop,
                              Simulation=type.simulations[6],
                              Vaccine="MMR2",
                              time=model.sol.0m.0r$model.sol[,"time"])

df.propor <- rbind(df.proport.0.V1,df.proport.0.V2,
                   df.proport.1.V1,df.proport.1.V2,
                   df.proport.2.V1,df.proport.2.V2,
                   df.proport.3.V1,df.proport.3.V2,
                   df.proport.4.V1,df.proport.4.V2,
                   df.proport.5.V1,df.proport.5.V2,
                   df.proport.6.V1,df.proport.6.V2)
df.propor$Simulation <- factor(df.propor$Simulation, levels = type.simulations)


# calling pdf printer
pdf(file = "measles_model_scenarios-3-vaccinated_proportions.pdf", width = 6, height = 4)

# visualizing results
df.propor %>%
  filter(time >= 18250) %>%
  mutate(time = time/365 + 1969) %>%
  ggplot() +
  geom_line(aes(x=time, y=Proportion, color=Simulation,
                linetype=Simulation), size = 1) +
  ggtitle("Proportion of protected individuals in the measles model") +
  ylab("Proportion")  + xlab("") +
  scale_x_continuous(breaks = seq(2020,2028, by = 2)) +
  theme(text = element_text(size=10)) +
  scale_color_manual(values = line.cols) +
  scale_linetype_manual(values = lines.typ) + 
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=3)) +
  facet_wrap(~Vaccine, scales = "free_y")

# turn off the pdf printer
dev.off()