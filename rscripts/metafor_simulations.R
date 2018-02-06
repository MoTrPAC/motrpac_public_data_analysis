# Meta analysis simulations
ncohorts = 10
nstudies = 5
yi = rnorm(ncohorts)
vi = rep(1,ncohorts)
S = rep(1:nstudies,(ncohorts/nstudies))
yi[S==1] = yi[S==1]+10
library(metafor)
rma.mv(yi,vi,random = ~1|S, mods = S)
