pedi<-c(0.425,0.35,0.38,0.86,0.67,0.771,0.78,0.59,0.61,0.66) # de 2009 en adelante fuente: https://ais.paho.org/imm/InfluenzaCoverageMap.asp
ges<-c(0,0,0,0,0.55,0.54,0.72,0.6,0.6,0.69) # 
adulm<-c(0.253,0.6,0.16,0.2,0.71,1,NA,0.5,NA,NA) #
year<-2009:2018
data1<-data.frame(year,pedi,ges,adulm)

# fitting models will be designed
start.year <- 2009
final.year <- 2018
number.years <- final.year - start.year + 1
data1$week <- seq(1, 52*number.years, 52)

par(mfrow=c(2,2))

# linear model for pediatric vaccine
pedi.model <- lm(formula = pedi ~ poly(week,7,raw = TRUE), data=data1)
coef(pedi.model)
t <- 1:(number.years*52)
pedi.fun <- (coef(pedi.model)[1] +  coef(pedi.model)[2]*t + coef(pedi.model)[3]*t^2 + coef(pedi.model)[4]*t^3 + coef(pedi.model)[5]*t^4 + coef(pedi.model)[6]*t^5 + coef(pedi.model)[7]*t^6 + coef(pedi.model)[8]*t^7)*(t<(9*52)) + mean(pedi[4:10])*(t>(9*52-1))# this is the model
plot(t/52+2009, pedi.fun, type = "l", main="Fitting pedriatic coverage",ylim=c(0,1),ylab="",xlab="")   
points(data1$week/52+2009, data1$pedi, col="red", pch=16)

# linear model for pregnant vaccine
ges.model <- mean(ges[5:10])
coef(pedi.model)
t <- 1:(number.years*52)
ges.fun <- ges.model*(t>(52*4))
plot(t/52+2009, ges.fun, type = "l", main="Fitting pregnant coverage",ylim=c(0,1),ylab="",xlab="")   
points(data1$week/52+2009, data1$ges, col="red", pch=16)

# linear model for ederly vaccine
adulm.model <- lm(formula = adulm ~ poly(week,4,raw = TRUE), data=data1)
coef(adulm.model)
t <- 1:(number.years*52)
adulm.fun <- (coef(adulm.model)[1] +  coef(adulm.model)[2]*t + coef(adulm.model)[3]*t^2 + coef(adulm.model)[4]*t^3 + coef(adulm.model)[5]*t^4) + mean(adulm[5],adulm[6],adulm[8])*(t>(7*52-1)) # this is the model
plot(t/52+2009, adulm.fun, type = "l", main="Fitting elderly coverage",ylim=c(0,1),ylab="",xlab="")    
points(data1$week/52+2009, data1$adulm, col="red", pch=16)

