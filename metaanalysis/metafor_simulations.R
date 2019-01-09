# Meta analysis simulations
ncohorts = 10
nstudies = 5
yi = rnorm(ncohorts)
vi = rep(1,ncohorts)
S = rep(1:nstudies,(ncohorts/nstudies))
yi[S==1] = yi[S==1]+10
library(metafor)
r1 = rma.mv(yi,vi,random = ~1|S, mods = S,intercept = F)
r2 = rma(yi,vi,mods = S)
r3 = rma.mv(yi,vi,random = ~1|S)
AIC(r2)
AIC(r3)
AIC(r1)
summary(r1)
confint(r1)
