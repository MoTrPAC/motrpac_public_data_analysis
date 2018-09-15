library(lme4)
library(simr)
library(longpower)

# Helper functions
# Simulate a linear mixed effects dataset
simulate_data<-function(n_t,sigma_between,sigma_within,effects_vec,n_subjects,effect_size,use_arima=F,arima_rho=0.5){
  intecepts = rnorm(n_subjects,sd = sigma_between)
  y = c()
  g = c()
  tp = c()
  for(j in 1:n_subjects){
    errs = rnorm(n_t,sd=sigma_within)
    if(use_arima){
      errs = arima.sim(list(order = c(1,0,0), ar = arima_rho), n = n_t,sd = sigma_within) 
    }
    effs = effect_size * effects_vec
    y = c(y,intecepts[j]+errs+effs)
    g = c(g,rep(j,n_t))
    tp = c(tp,0:(n_t-1))
  }
  d = data.frame(y,g,t=tp)
  return(d)
}
# Plot the trajectories of the subjects
plot_longi<-function(d){
  plot(y=d$y,d$t,type="n",ylab="abundance",xlab="time")
  points(y=d$y,d$t,type="p",pch=3)
  for(i in unique(d$g)){lines(y=d$y[d$g==i],d$t[d$g==i],type="l",lty=2)} 
  lines(lowess(y=d$y,d$t),lwd=3,col="yellow") 
}

# Example dataset:
n_t = 5 # one for pre and then four time points
sigma_between = 2 # random effect standard deviation
n_subjects = 20
effects_vec = c(0,0.25,1,0.25,0)
sigma_within = 1
effect_size = 1
d = simulate_data(n_t,sigma_between,sigma_within ,effects_vec,n_subjects,effect_size)
plot_longi(d)

n_t = 5 # one for pre and then four time points
sigma_between = 1 # random effect standard deviation
n_subjects = 100
effects_vec = c(0,0.25,1,0.25,0)
sigma_within = 1
effect_size = 0.25
nsim=1000

d = simulate_data(n_t,sigma_between,sigma_within ,effects_vec,n_subjects,effect_size)
plot_longi(d)

library(lme4)
model = lmer(y~factor(t)+(1|g),data=d)
model = lmer(y~ordered(t)+(1|g),data=d)
summary(model)

# use simr
simr_model = model
# specify desired effect sizes
for(j in 2:n_t){
fixef(simr_model)[paste("factor(t)",j-1,sep="")] = effects_vec[j]*effect_size
}
# Analysis at a range of sample sizes
model3 <- extend(simr_model, along="g", n=600)
pc3 <- powerCurve(model3, along="g",test=fixed(paste("factor(t)",2,sep=""),"z"),alpha=0.001,nsim=50)
plot(pc3,xlab="Number of subjects")
powerSim(simr_model,test=fixed(paste("factor(t)",1,sep=""),"z"),alpha=0.001,nsim=50)





