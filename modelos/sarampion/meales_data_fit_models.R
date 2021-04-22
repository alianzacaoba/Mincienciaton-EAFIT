#############################
# libraries
#############################

library(dplyr)


#############################
# loading measles data
#############################

data.ori <- data.table::fread(input = "population-birth_death_rates-vaccov-1960-2019.csv",
                          header=TRUE, sep=",", na.strings = "",
                          col.names = c("year","population","birth.rate","death.rate",
                                        "vaccov1","vaccov2","cases"))
data.ori$death.rate <- data.ori$death.rate*1e-3 # annual & crude (per 1000 people)
data.ori$birth.rate <- data.ori$birth.rate*1e-3 # annual & crude (per 1000 people)
data <- filter(data.ori, year >= 1970 & year <= 2018)


#############################
# defining temporal frame
#############################

# available temporal frame based on which
# fitting models will be designed
start.year <- 1970
final.year <- 2018
number.years <- final.year - start.year + 1
data$days <- seq(1, 365*number.years, 365)


#############################
# fitting models
#############################

# linear model for population
pop.model <- lm(formula = population ~ days, data=data)
coef(pop.model)

# linear model for birth rate
birth.model <- lm(formula = birth.rate ~ days, data=data)
coef(birth.model)

# linear model for death rate
death.model <- lm(formula = death.rate ~ poly(days,7,raw = TRUE), data=data)
coef(death.model)


# bounded exponential model for vaccine 1
vacc1.intro.year <- 1980
data.vacc1 <- filter(data, year>=vacc1.intro.year)
vacc1.model <- nls(vaccov1 ~ 0.95*(1-exp(-k*(days-3284)/365)),
                   data = data.vacc1, 
                   algorithm="port",
                   start=c(k=0.1),lower=c(k=0.07), upper=c(k=0.2))
coef(vacc1.model)


# bounded exponential model for vaccine 2
vacc2.intro.year <- 1997
data.vacc2 <-  filter(data, year>=vacc2.intro.year)
vacc2.model <- nls(vaccov2 ~ (0.85-exp(-k*(days-9489)/365)),
                   data = data.vacc2, 
                   algorithm="port",
                   start=c(k=0.35),lower=c(k=0.2), upper=c(k=0.4))
coef(vacc2.model)


#############################
# visualizing results
#############################

# vector of days
# (from 1st day of 1980 to the last day of 2018)
# recall: number.years <- final.year - start.year + 1
t <- 1:(number.years*365)


# to plot all the graphics in the same figure
layout(matrix(1:6, 2, 3, byrow = FALSE))


# plotting population data (red points) and the fitted model (solid line)
pop.fun <- t*coef(pop.model)[2] + coef(pop.model)[1] # this is the model
plot(t, pop.fun, type = "l", main="Fitting population data")  
points(data$days, data$population, col="red", pch=16)


# plotting birth rate data (red points) and the fitted model (solid line)
birth.fun <- t*coef(birth.model)[2] + coef(birth.model)[1] # this is the model
plot(t, birth.fun, type = "l", main="Fitting birth rate data")  
points(data$days, data$birth.rate, col="red", pch=16)


# plotting death rate data (red points) and the fitted model (solid line)
death.fun <- coef(death.model)[1] +  coef(death.model)[2]*t + coef(death.model)[3]*t^2 + coef(death.model)[4]*t^3 + coef(death.model)[5]*t^4 + coef(death.model)[6]*t^5 + coef(death.model)[7]*t^6 + coef(death.model)[8]*t^7 # this is the model
plot(t, death.fun, type = "l", main="Fitting death rate data")  
points(data$days, data$death.rate, col="red", pch=16)


# plotting 1st dose coverage data (red points) and the fitted model (solid line)
vacc1.fun <- (t < 3651)*0 + 
             (t >= 3651)*(0.95*(1-exp(-coef(vacc1.model)*(t-3284)/365))) # this is the model
plot(t, vacc1.fun, type = "l", main="Fitting 1st dose coverage data")  
points(data.vacc1$days, data.vacc1$vaccov1, col="red", pch=16)


# plotting 2nd dose coverage data (red points) and the fitted model (solid line)
vacc2.fun <- (t < 9856)*0 + 
             (t >= 9856)*(0.85-exp(-coef(vacc2.model)*(t-9489)/365)) # this is the model
plot(t, vacc2.fun, type = "l", main="Fitting 2nd dose coverage data")  
points(data.vacc2$days, data.vacc2$vaccov2, col="red", pch=16)


#########################################################
# report figure
#########################################################

# calling pdf printer
pdf(file = "common_fitted_models.pdf", width = 6, height = 7)

# to plot all the graphics in the same figure
layout(matrix(1:6, 3, 2, byrow = TRUE))

# plotting population data (red points) and the fitted model (solid line)
pop.fun <- t*coef(pop.model)[2] + coef(pop.model)[1] # this is the model
plot(t, pop.fun/1e6, type = "l", xlab = "Days", ylab = "Population (millions)",
     main="Population time series and fitted model",
     cex=0.5, cex.main=0.9)  
points(data$days, data$population/1e6, pch=16,
       col=rgb(red=1, green=0, blue=0, alpha=0.5))

# plotting birth rate data (red points) and the fitted model (solid line)
birth.fun <- t*coef(birth.model)[2] + coef(birth.model)[1] # this is the model
plot(t, birth.fun, type = "l", xlab = "Days", ylab = "Birth rate",
     main="Birth rate time series and fitted model",
     cex=0.5, cex.main=0.9)   
points(data$days, data$birth.rate, pch=16,
       col=rgb(red=1, green=0, blue=0, alpha=0.5))

# plotting death rate data (red points) and the fitted model (solid line)
death.fun <- coef(death.model)[1] +  coef(death.model)[2]*t + coef(death.model)[3]*t^2 + coef(death.model)[4]*t^3 + coef(death.model)[5]*t^4 + coef(death.model)[6]*t^5 + coef(death.model)[7]*t^6 + coef(death.model)[8]*t^7 # this is the model
plot(t, death.fun, type = "l", xlab = "Days", ylab = "Death rate",
     main="Death rate time series and fitted model",
     cex=0.5, cex.main=0.9)    
points(data$days, data$death.rate, pch=16,
       col=rgb(red=1, green=0, blue=0, alpha=0.5))

# plotting 1st dose coverage data (red points) and the fitted model (solid line)
vacc1.fun <- (t < 3651)*0 + 
  (t >= 3651)*(0.95*(1-exp(-coef(vacc1.model)*(t-3284)/365))) # this is the model
plot(t, vacc1.fun, type = "l", xlab = "Days", ylab = "MMR1 coverage",
     main="MMR1 coverage time series and fitted model",
     cex=0.5, cex.main=0.9)    
points(data.vacc1$days, data.vacc1$vaccov1, pch=16,
       col=rgb(red=1, green=0, blue=0, alpha=0.5))

# plotting 2nd dose coverage data (red points) and the fitted model (solid line)
vacc2.fun <- (t < 9856)*0 + 
  (t >= 9856)*(0.85-exp(-coef(vacc2.model)*(t-9489)/365)) # this is the model
plot(t, vacc2.fun, type = "l", xlab = "Days", ylab = "MMR2 coverage",
     main="MMR2 coverage time series and fitted model",
     cex=0.5, cex.main=0.9)   
points(data.vacc2$days, data.vacc2$vaccov2, pch=16,
       col=rgb(red=1, green=0, blue=0, alpha=0.5))

# turn off the pdf printer
dev.off()