### THE FOLLOWING COMMAND SOURCES THE R CODE WITH THE SAMPLER IMPLEMENTATIONS ###
source("https://raw.githubusercontent.com/gzanella/TGS/master/functions_for_BVS.R")

###### SPECIFYING THE MODEL AND SIMULATING DATA #######
####
n<-100 # n.observations
p<-100 # n.covariates
c<-10^3 #prior hyperparameter for the covariance matrix
prior_p_incl<-5/p #prior prob of inclusion
hyper_par<-simulate_data(n=n,p=p,c=c,SNR=3,scenario=1) # creates list to be passed to samplers

###### RUNNING THE SAMPLERS #######
####
T<-50000 #n.iterations
burn_in<-5000
thin_GS<-2
output_GS<-GS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in,thin = thin_GS)
output_TGS<-TGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in)
output_wTGS<-wTGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in)


##### PLOTTING THE OUTPUT ##########
####

par(mfrow=c(1,2))
col_GS<-"black";col_TGS<-"blue";col_wTGS<-"red"
# PLOT RUNNING ESTIMATES OF Pr(gamma_1=1|Y) FOR DIFFERENT SAMPLERS
plot(output_TGS$est_inclusion_probs_1/cumsum(output_TGS$sample_weights),
     ylim=c(0,1),type="l",col=col_TGS,xlab="Sampler iteration", ylab="Posterior inclusion prob.")
lines(cumsum(output_GS$gamma_1[(burn_in+1):(burn_in+T)])/c(1:T),
      col=col_GS)
lines(output_wTGS$est_inclusion_probs_1/cumsum(output_wTGS$sample_weights),type="l",col=col_wTGS)
title("Variable 1")
# PLOT RUNNING ESTIMATES OF Pr(gamma_2=1|Y) FOR DIFFERENT SAMPLERS
plot(output_TGS$est_inclusion_probs_2/cumsum(output_TGS$sample_weights),
     ylim=c(0,1),type="l",col=col_TGS,xlab="Sampler iteration", ylab="Posterior inclusion prob.")
lines(cumsum(output_GS$gamma_2[(burn_in+1):(burn_in+T)])/c(1:T),col=col_GS)
lines(output_wTGS$est_inclusion_probs_2/cumsum(output_wTGS$sample_weights),type="l",col=col_wTGS)
title("Variable 2")
legend(x='topright',legend=c('GS','TGS','wTGS'),lty=c(1,1,1),col=c(col_GS,col_TGS,col_wTGS))
par(mfrow=c(1,2))