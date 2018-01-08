library(data.table)

simLogData <- function(x,a=1,sd=3,samples=10){
  foo=list()
  for(z in x){
    foo[[z]]=data.table(x=z,y=rnorm(samples,mean = a*log(z),sd = sd))
  }
  return(rbindlist(foo))
}

## Simulate data 
dat=simLogData(1:10,a=5,sd=1)

## Plot data
plot(dat$x,dat$y)
## Fit log model
log.mod=glm(formula = y~log(x),family = "gaussian", data = dat)
## Print summary of model
summary(log.mod)
## Predict values from model for plotting
pred.tab=data.table(x=seq(min(dat$x),max(dat$x),0.1))
pred.tab[,y:=predict(object = log.mod,newdata = pred.tab)]
lines(x=pred.tab$x,y = pred.tab$y)

## Fit linear model
lin.mod=glm(formula = y~x,family = "gaussian", data = dat)
## Print summary of model
summary(lin.mod)

## Compute difference in AIC (smaller is better)
## Negative means the the log model is better
## The models are usually different if AIC >= 2
log.mod$aic - lin.mod$aic 

## Also if you want to calculate the velocity at a point just use this function
velocity <- function(x,mod){
  coef(mod)[2]/x
}

velocity(1,log.mod)
velocity(1,lin.mod)
