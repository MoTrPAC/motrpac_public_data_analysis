library(metafor)
# Meta analysis simulations: check different time models
ncohorts = 10
nstudies = 5
ntimes = 5
N = ntimes*ncohorts
yi = rnorm(N)
vi = rep(1,N)
times = rep(1:ntimes,ncohorts)
cohorts = sort(rep(1:ncohorts,ntimes))
# add true time-dep effect
yi[times<3] = yi[times<3]+5
REs = rnorm(ncohorts)
for(j in 1:ncohorts){
  curr_re = rnorm(1)
  yi[cohorts==j] = yi[cohorts==j]+curr_re
}
times = ordered(times)
d=data.frame(yi,vi,times=times,cohorts)
r1 = rma.mv(yi,vi,random = ~1|times,data=d)
fs = fitstats(r1)
get_df_aic(fs[1,1],fs[3,1])

r1 = rma.mv(yi,vi,random = ~1|cohorts,mods=~times,data=d)
fs = fitstats(r1)
fs
get_df_aic(fs[1,1],fs[3,1])

r1 = rma.mv(yi,vi,random = ~times|cohorts,mods=~times,data=d)
fs = fitstats(r1)
fs
get_df_aic(fs[1,1],fs[3,1])

r1 = rma(yi,vi,mods=~times,data=d)
fs = fitstats(r1)
get_df_aic(fs[1,1],fs[3,1])

r1 = rma(yi,vi,data=d)
fs = fitstats(r1)
get_df_aic(fs[1,1],fs[3,1])

get_df_aic<-function(loglik,aic){
  k = (aic+2*loglik)/2
  return(k)
}
get_df_bic<-function(loglik,bic,n){
  lnk = (bic+2*loglik)
  return(lnkk/log(n))
}



