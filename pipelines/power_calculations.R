library(lme4)
library(simr)
library(longpower)

# Helper functions
# Simulate a linear mixed effects dataset
simulate_data<-function(n_t,sigma_between,sigma_within,effects_vec,
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
# Plot the trajectories of the subjects
plot_longi<-function(d){
  plot(y=d$y,d$t,type="n",ylab="abundance",xlab="time")
  points(y=d$y,d$t,type="p",pch=3)
  for(i in unique(d$g)){lines(y=d$y[d$g==i],d$t[d$g==i],type="l",lty=2)} 
  lines(lowess(y=d$y,d$t),lwd=3,col="yellow") 
}

# Example dataset:
n_t = 5 # one for pre and then four time points
sigma_between = 1 # random effect standard deviation
n_subjects = 100
effects_vec = c(0,0.25,1,0.25,0)
sigma_within = 0.1
effect_size = 1
d = simulate_data(n_t,sigma_between,sigma_within ,effects_vec,n_subjects,effect_size)
plot_longi(d)

library(lme4);library(simr);library(gplots)
#' Power analysis that treats time as simple factors
#' 
#'  @param d A data frame of simulated data. For example, output of simulate_data
#'  @param effects_var A numeric vector specifying relative weights for a trajectory
#'  @param effect_size A number. the effect size to multiply effects_var by
#'  @param max_m A number.  The maximal numer of subjects to consider in calculations
#'  @param nsim A number. Number of simulations for simr
#'  @param alpha A number. The significance level for the power calculations
#'  @param tp A number. The time point of interest. The power calculation focuses on 
#'  the ability to detect the effect in this time point
#'  
#'  @return A power curve
get_simple_power_plot<-function(d,effects_vec,effect_size,
                                max_n=700,nsim=10,alpha=0.001,tp=1,...){
  simr_model = lmer(y~factor(t)+(1|g),data=d)
  # specify desired effect sizes
  for(j in 2:n_t){
    fixef(simr_model)[paste("factor(t)",j-1,sep="")] = effects_vec[j]*effect_size
  }
  # Analysis at a range of sample sizes
  model3 <- extend(simr_model, along="g", n=max_n)
  pc3 <- powerCurve(model3, along="g",
                    test=fixed(paste("factor(t)",tp,sep=""),"z"),
                    alpha=alpha,nsim=nsim,...)
  # plot(pc3,xlab="Number of subjects")
  return(pc3)
}

#' Power analysis that uses polynomials to model time trajectories.
#' 
#'  @param d A data frame of simulated data. For example, output of simulate_data
#'  @param effects_var A numeric vector specifying relative weights for a trajectory
#'  @param effect_size A number. the effect size to multiply effects_var by
#'  @param max_m A number.  The maximal numer of subjects to consider in calculations
#'  @param nsim A number. Number of simulations for simr
#'  @param alpha A number. The significance level for the power calculations
#'  @param poly_deg A number. The polynomial degree. The power calculation focuses on 
#'  the ability to detect the effect in this trajectory type
#'  
#'  @return A power curve
get_poly_power_plot<-function(d,effects_vec,effect_size,
                              max_n=700,nsim=10,alpha=0.001,poly_deg=2,...){
  n_t = length(unique(d$t)) 
  pp = poly(0:(n_t-1),2)
  yy = effects_vec*effect_size
  new_effects_vec = lm(yy~pp)$coefficients[-1]
  rownames(pp) = as.character(0:(n_t-1))
  d_pp = pp[as.character(d$t),]
  rownames(d_pp) = NULL
  d = data.frame(d,d_pp)
  model = lmer(y~X1+X2+(1|g),data=d)
  fnames = summary(model)$coefficients
  fnames = rownames(fnames)[-1]
  simr_model = model
  # specify desired effect sizes
  for(j in 1:length(fnames)){fixef(simr_model)[fnames[j]] = new_effects_vec[j]}
  # Analysis at a range of sample sizes
  model3 <- extend(simr_model, along="g", n=max_n)
  pc3 <- powerCurve(model3, along="g",
                    test=fixed(fnames[poly_deg],"z"),
                    alpha=alpha,nsim=nsim,...)
  # plot(pc3,xlab="Number of subjects")
  return(pc3)
}

n_t = 5 # one for pre and then four time points
sigma_between = 1 # random effect standard deviation
n_subjects = 100
effects_vec = c(0,0.25,1,0.25,0)
sigma_within = 0.5
effect_size = 0.25
nsim=100
d = simulate_data(n_t,sigma_between,sigma_within ,effects_vec,n_subjects,effect_size)
par(mfrow=c(2,2))
ps1 = get_simple_power_plot(d,effects_vec,effect_size,max_n=700,nsim=10,alpha=0.001,tp=2)
ps2 = get_simple_power_plot(d,effects_vec,effect_size,max_n=700,nsim=10,alpha=0.001,tp=1)
ps3 = get_poly_power_plot(d,effects_vec,effect_size,max_n=700,nsim=10,alpha=0.001)
plot_ci_results(list("S:tp=1"=ps1,"S:tp=2"=ps2,"P:q"=ps3))

plot(ps1)
lines(ps2)

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

n_subjects = 120
effects_vec = c(0,0.25,1,0.25,0)
sigma_within_range = c(0.25,0.5,1)
effect_size_range = seq(0.1,1,by=0.1)
nsim=10
res_tp2 = list();res_q = list()
for(sigma_within in sigma_within_range){
  m1 = c();m2=c()
  for (effect_size in effect_size_range){
    d = simulate_data(n_t,sigma_between,sigma_within ,effects_vec,n_subjects,effect_size)
    ps_tp = get_simple_power_plot(d,effects_vec,effect_size,max_n=n_subjects,
                                nsim=nsim,alpha=0.001,tp=2,breaks=n_subjects)
    ps_tp = summary(ps_tp)[1,]
    ps_tp[1] = effect_size
    m1 = rbind(m1,ps_tp)
    ps_q = get_poly_power_plot(d,effects_vec,effect_size,max_n=n_subjects,
                              nsim=nsim,alpha=0.001,breaks=n_subjects)
    ps_q = summary(ps_q)[1,]
    ps_q[1] = effect_size
    m2 = rbind(m2,ps_q)
  }
  res_tp2[[as.character(sigma_within)]] = m1
  res_q[[as.character(sigma_within)]] = m2
}
plot_ci_results(res_tp2,xlab="Peak effect size")
