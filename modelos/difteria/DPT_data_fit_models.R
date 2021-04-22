#############################
# libraries
#############################

library(dplyr)


#############################
# loading measles data
#############################

data <- data.table::fread(input = "dpt3_vaccine_coverage_OMS.csv",
                          header=TRUE, sep=";", na.strings = "")



#############################
# defining temporal frame
#############################

# available temporal frame based on which
# fitting models will be designed
start.year <- 1970
final.year <- 2019
number.years <- final.year - start.year + 1
data$days <- seq(1, 365*number.years, 365)


#############################
# fitting models
#############################

# bounded exponential model for vaccine 1
vacc1.intro.year <- 1980
data.vacc1 <- filter(data, Year >= vacc1.intro.year)
vacc1.model <- nls(DPT3_percentage ~ 0.92*(1-exp(-k*(days)/365)),
                   data = data.vacc1, 
                   algorithm="port",
                   start=c(k=0.1),lower=c(k=0.07), upper=c(k=0.2))
coef(vacc1.model)



vc.model <- nls(formula = DPT3_percentage ~ 1-exp(-k*(days)/365), data=datos,
                   algorithm="port",
                   start=c(k=0.1),lower=c(k=0.07), upper=c(k=0.2))
coef(vc.model)

#############################
# visualizing results
#############################

# vector of days
# (from 1st day of 1980 to the last day of 2018)
# recall: number.years <- final.year - start.year + 1
t <- 1:(number.years*365)

# plotting 1st dose coverage data (red points) and the fitted model (solid line)
vacc1.fun <- (0.92*(1-exp(-coef(vacc1.model)*(t)/365))) # this is the model
plot(t, vacc1.fun, type = "l", main="Fitting DTP3 coverage data")  
points(data.vacc1$days, data.vacc1$DPT3_percentage, col="red", pch=16)

# plotting death rate data (red points) and the fitted model (solid line)
vc.fun <- 0.92*(1-exp(-coef(vc.model)*(t)/365))  # this is the model
plot(t, vc.fun, type = "l", main="Fitting vaccine coverage data")  
points(datos$days, datos$DPT3_percentage, col="red", pch=16)
