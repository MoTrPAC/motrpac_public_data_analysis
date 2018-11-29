library(lme4)
library(simr)
library(longpower)

# Helper functions
# Simulate a linear mixed effects dataset
# Focus here is on creating a dataset with a temporal component
# and subject-specific random effects.
simulate_initial_dataset<-function(n_t,sigma_between,sigma_within,effects_vec,
                        n_subjects,effect_size,use_arima=F,arima_rho=0.5){
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

# Here we create subject specific covariates such as age and group.
# cov is a vector with a value for each subject.
# For each subject an effect is added based on the covariate,
# but not in time point zero in order to model differential
# response and not just an overall fixed effect (that will keep the same response)
add_covariate_to_data<-function(d,cov,effect=1){
  v = rep(1,nrow(d))
  new_y = d$y
  for(i in 1:length(cov)){
    curr_inds = d$g==i
    v[curr_inds] = cov[i]
    curr_inds = curr_inds  & d$t>0
    new_y[curr_inds] = new_y[curr_inds] + cov[i]*effect
  }
  d = data.frame(y=new_y,g=d$g,t=d$t,v=v)
  return(d)
}

# Plot the trajectories of the subjects
plot_longi<-function(d,...){
  plot(y=d$y,d$t,type="n",xlab="time",...)
  points(y=d$y,d$t,type="p",pch=3)
  for(i in unique(d$g)){lines(y=d$y[d$g==i],d$t[d$g==i],type="l",lty=2)} 
  lines(lowess(y=d$y,d$t),lwd=3,col="yellow") 
}

library(lme4);library(simr);library(gplots)
#' Power analysis that treats time as simple factors
#' 
#'  @param d A data frame of simulated data. For example, output of simulate_data
#'  @param effects_var A numeric vector specifying relative weights for a trajectory
#'  @param effect_size1 A number. the effect size to multiply effects_var by
#'  @param effect_size2 The effect size of the additional covariate
#'  @param max_m A number.  The maximal numer of subjects to consider in calculations
#'  @param nsim A number. Number of simulations for simr
#'  @param alpha A number. The significance level for the power calculations
#'  
#'  @return A power curve
get_simple_power_plot<-function(d,effects_vec,effect_size1,n_t,effect_size2=1,
                                max_n=300,nsim=10,alpha=0.001,...){
  simr_model = lmer(y~factor(t) + v +(1|g),data=d)
  # specify desired effect sizes
  for(j in 2:n_t){
    fixef(simr_model)[paste("factor(t)",j-1,sep="")] = effects_vec[j]*effect_size1
  }
  fixef(simr_model)["v"] = effect_size2
  # Analysis at a range of sample sizes
  model3 <- extend(simr_model, along="g", n=max_n)
  pc1 <- powerCurve(model3, along="g",
                    test=fixed("factor(t)1","z"),alpha=alpha,nsim=nsim)
  pc2 <- powerCurve(model3, along="g",
                    test=fixed("v","z"),alpha=alpha,nsim=nsim)
  # plot(pc3,xlab="Number of subjects")
  return(list(pc1,pc2))
}

library(gplots)
#' Auxiliary function for plotting power curves with errors.
#' 
#'  @param l A list of powerCurve objects
#'  @param cols A vector of colors, length(cols)>=length(l)
#'  @param cols A vector of pch codes, length(pchs)>=length(l)
#'  @param ... Additional paramaters for legend
plot_ci_results<-function(l,cols = c("red","green","blue"),pchs=20:24,
                          xlab="Number of subjects",...){
  if(class(l[[1]])[1]=="powerCurve"){l = lapply(l,summary)}
  plot(l[[1]][,1], l[[1]][,4],ylim=c(0,1.2),type="b",col=cols[1],pch=pchs[1],
       xlab=xlab,ylab="power")
  for(j in 1:length(l)){
    x = l[[j]][,1]
    ui = pmin(1,l[[j]][,6])
    li = pmax(0,l[[j]][,5])
    arrows(x, li, x, ui, length=0.05, angle=90, code=3,col=cols[j])  
    if(j>1){
      lines(l[[j]][,1], l[[j]][,4],type="b",col=cols[j],pch=pchs[j])
    }
  }
  abline(h=0.8,lty=2)
  legend(x="topleft",legend=names(l),pch=pchs,col=cols,lwd=2,ncol=length(cols),...)
}

# Example dataset:
n_t = 2 # one for pre and then four time points
sigma_between = 1 # random effect standard deviation
sigma_within = 1
n_subjects = 50
effects_vec = c(0,1)

# Test and plot
d = simulate_initial_dataset(n_t,sigma_between,sigma_within ,effects_vec,n_subjects,effect_size)
plot_longi(d)
age_cov = sort(rep(0:1,n_subjects/2))
newd = add_covariate_to_data(d,age_cov)
curr_sim = get_simple_power_plot(newd,effects_vec,effect_size1=effect_size,n_t=2,
                                 effect_size2=effect_size,max_n=150,nsim=50,alpha=0.001)
par(mfrow=c(2,2))
plot_longi(newd[newd$v==0,],main="Older",xaxt="n",ylab="Max Mito Capacity")
axis(1, xaxp=c(0, 1, 1), las=1)
plot_longi(newd[newd$v==1,],main="Middle-aged",xaxt="n",ylab="Max Mito Capacity")
axis(1, xaxp=c(0, 1, 1), las=1)
plot(curr_sim[[1]],xlab="Number of subjects")
title("General trend")
plot(curr_sim[[2]],xlab="Number of subjects")
title("Age-specific effect")

sigma_within_range = c(0.25,0.5,1)
effect_size_range = seq(0.1,1,by=0.1)
nsim=10
sim_results_main = list()
sim_results_age = list()
for(sigma_within in sigma_within_range){
  m1 = c();m2=c()
  for (effect_size in effect_size_range){
    d = simulate_initial_dataset(n_t,sigma_between,sigma_within ,effects_vec,n_subjects,effect_size)
    age_cov = sort(rep(0:1,n_subjects/2))
    newd = add_covariate_to_data(d,age_cov,effect_size)
    
    curr_sim = get_simple_power_plot(newd,effects_vec,effect_size1=effect_size,n_t=2,
                                     effect_size2=effect_size,max_n=n_subjects,nsim=nsim,alpha=0.001)
    ps_tp = summary(curr_sim[[1]])[1,]
    ps_tp[1] = effect_size
    m1 = rbind(m1,ps_tp)
    ps_tp = summary(curr_sim[[2]])[1,]
    ps_tp[1] = effect_size
    m2 = rbind(m2,ps_tp)
  }
  sim_results_main[[as.character(sigma_within)]] = m1
  sim_results_age[[as.character(sigma_within)]] = m2
}
plot_ci_results(sim_results_main,xlab="Peak effect size")
