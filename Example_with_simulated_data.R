setwd('C:/Users/ZanellaG/Bitbucket/TGS/')
source('functions_for_BVS.R')

###### SPECIFYING THE MODEL AND SIMULATING DATA #######
####
n<-100
p<-100
c<-10^3
prior_p_incl<-5/p
hyper_par<-simulate_data(n=n,p=p,c=c,SNR=3,scenario=1)

###### RUNNING THE SAMPLERS #######
####
T<-50000
burn_in<-5000
output_TGS<-TGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in)
thin_GS<-5
output_GS<-GS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in,thin = thin_GS)
output_wTGS<-wTGS(p=p,hyper_par=hyper_par,T=T,burn_in=burn_in)


##### PLOTTING THE OUTPUT ##########
####

par(mfrow=c(1,2))
col_GS<-"black";col_TGS<-"blue";col_wTGS<-"red"
# plot running estimates of Pr(gamma_1=1|Y) for different samplers
plot(output_TGS$est_inclusion_probs_1/cumsum(output_TGS$sample_weights),
     ylim=c(0,1),type="l",col=col_TGS,xlab="Sampler iteration", ylab="Posterior inclusion prob.")
lines(cumsum(output_GS$gamma_1[(burn_in+1):(burn_in+T)])/c(1:T),
      col=col_GS)
lines(output_wTGS$est_inclusion_probs_1/cumsum(output_wTGS$sample_weights),type="l",col=col_wTGS)
title("Variable 1")
# plot running estimates of Pr(gamma_2=1|Y) for different samplers
plot(output_TGS$est_inclusion_probs_2/cumsum(output_TGS$sample_weights),
     ylim=c(0,1),type="l",col=col_TGS,xlab="Sampler iteration", ylab="Posterior inclusion prob.")
lines(cumsum(output_GS$gamma_2[(burn_in+1):(burn_in+T)])/c(1:T),col=col_GS)
lines(output_wTGS$est_inclusion_probs_2/cumsum(output_wTGS$sample_weights),type="l",col=col_wTGS)
title("Variable 2")
legend(x='topright',legend=c('GS','TGS','wTGS'),lty=c(1,1,1),col=c(col_GS,col_TGS,col_wTGS))
par(mfrow=c(1,2))