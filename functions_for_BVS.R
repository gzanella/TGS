#### MAIN FUNCTIONS FOR SAMPLERS ######
###
GS<-function(p,gamma_start=NULL,
             T,burn_in=0,thin=1,
             hyper_par=NULL,vars_selected=c(1,2)){ 
  ## random scan Gibbs Sampler for Bayesian variable selection problems
  
  
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  indices_sequence<-rep(NA,burn_in+T)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    ## perform random scan Gibbs Sampling step
    for (iter in 1:thin){
      j<-sample.int(n = p,size = 1)
      fc_j<-single_full_cond(j=j,gamma=gamma,hyper_par=hyper_par)$fc_j
      # ap<-1-fc_j## not metropolized
      ap<-(1-fc_j)/fc_j## metropolized
      if(runif(1) < ap){
        gamma[j]<-1-gamma[j]
        indices_sequence[t]<-j
      }
    }
    ## update inclusion estimators
    if(t>burn_in){
      est_inclusion_probs<-est_inclusion_probs+gamma/T
    }
    ## store indices sequence for output analysis
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
  }
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma))
}

TGS<-function(p,gamma_start=NULL,
              T,burn_in=0,
              hyper_par=NULL,vars_selected=c(1,2)){
  ## TGS algorithm for Bayesian variable selection problems
  
  
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  est_inclusion_probs_1<-rep(NA,T)
  est_inclusion_probs_2<-rep(NA,T)
  indices_sequence<-rep(NA,burn_in+T)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  sample_weights<-rep(NA,T)
  ## compute full conditionals
  output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=NULL)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    ## perform Tempered Gibbs Sampling step
    if(length(which(output_fc$fc==0))>0){
      j<-sample(which(output_fc$fc==0),size = 1)
    }else{
      j<-which(cumsum(1/output_fc$fc)>(runif(1)*sum(1/output_fc$fc)))[1]
    }
    gamma[j]<-1-gamma[j]
    # gamma[j]<-rbinom(1,1,0.5)
    ## compute full conditionals and update Rao-Blackwelied estimator
    output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=output_fc$stored)
    if(t>burn_in){
      sample_weights[t-burn_in]<-(p/2)/sum(1/output_fc$fc)
      est_inclusion_probs<-est_inclusion_probs+sample_weights[t-burn_in]*
        (gamma*output_fc$fc+(1-gamma)*(1-output_fc$fc))
      est_inclusion_probs_1[t-burn_in]<-est_inclusion_probs[vars_selected[1]]
      est_inclusion_probs_2[t-burn_in]<-est_inclusion_probs[vars_selected[2]]
    }
    ## store indices sequence for output analysis
    indices_sequence[t]<-j
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
  }
  est_inclusion_probs<-est_inclusion_probs/sum(sample_weights)
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              sample_weights=sample_weights,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma,
              est_inclusion_probs_1=est_inclusion_probs_1,
              est_inclusion_probs_2=est_inclusion_probs_2))
}

wTGS<-function(p,gamma_start=NULL,
              T,burn_in=0,
              hyper_par=NULL,vars_selected=c(1,2)){
  ## wTGS algorithm for Bayesian variable selection problems
  
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  est_inclusion_probs_1<-rep(NA,T)
  est_inclusion_probs_2<-rep(NA,T)
  indices_sequence<-rep(NA,burn_in+T)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  sample_weights<-rep(NA,T)
  ## compute full conditionals
  output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=NULL)
  flip_rates<-compute_flip_rates(output_fc=output_fc,gamma=gamma)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
   ## perform Markov chain transition step
   if(length(which(output_fc$fc==0))>0){
      j<-sample(which(output_fc$fc==0),size = 1)
    }else{
      j<-which(cumsum(flip_rates)>(runif(1)*sum(flip_rates)))[1]
    }
    gamma[j]<-1-gamma[j]
    ## compute full conditionals and update Rao-Blackwelied estimator
    output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=output_fc$stored)
    flip_rates<-compute_flip_rates(output_fc=output_fc,gamma=gamma)
    # flip_rates[which(gamma==0)]<-(1-output_fc$fc[which(gamma==0)])/output_fc$fc[which(gamma==0)]
    # flip_rates[which(gamma==1)]<-1
    if(t>burn_in){
      sample_weights[t-burn_in]<-(p/2)/sum(flip_rates)
      est_inclusion_probs<-est_inclusion_probs+sample_weights[t-burn_in]*
        (gamma*output_fc$fc+(1-gamma)*(1-output_fc$fc))
      est_inclusion_probs_1[t-burn_in]<-est_inclusion_probs[vars_selected[1]]
      est_inclusion_probs_2[t-burn_in]<-est_inclusion_probs[vars_selected[2]]
    }
    ## store indices sequence for output analysis
    indices_sequence[t]<-j
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
  }
  est_inclusion_probs<-est_inclusion_probs/sum(sample_weights)
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              sample_weights=sample_weights,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma,
              est_inclusion_probs_1=est_inclusion_probs_1,
              est_inclusion_probs_2=est_inclusion_probs_2))
}


GS_RB<-function(p,gamma_start=NULL,
               T,burn_in=0,thin=1,
               hyper_par=NULL,vars_selected=c(1,2)){
  ## deterministic scan Gibbs Sampler with Rao-Blackwellization for Bayesian variable selection problems
  
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    for (j in 1:p){
      fc_j<-single_full_cond(j=j,gamma=gamma,hyper_par=hyper_par)$fc_j
      # ap<-1-fc_j## not metropolized
      ap<-(1-fc_j)/fc_j## metropolized
      ## update inclusion estimators
      if(t>burn_in){
        est_inclusion_probs[j]<-est_inclusion_probs[j]+
          (gamma[j]*fc_j+(1-gamma[j])*(1-fc_j))/T
      }
      if(runif(1) < ap){
        gamma[j]<-1-gamma[j]
      }
    }
    ## store indices sequence for output analysis
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
  }
  return(list(est_inclusion_probs=est_inclusion_probs,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma))
}

HBS<-function(p,gamma_start=NULL,full_cond,
              T,burn_in=0,
              hyper_par=NULL){
  ## Hamming Ball Sampler for Bayesian variable selection problems

  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  indices_sequence<-matrix(rep(NA,2*(burn_in+T)),nrow=2)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    ## perform Hamming ball sampler move (i.e. sample U|X and X|U)
    j_aux<-sample.int(p,1)
    gamma[j_aux]<-1-gamma[j_aux]
    output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=NULL)
    probs<-(1-output_fc$fc)/output_fc$fc
    if(length(which(probs==Inf))>0){
      j<-sample(which(probs==Inf),size = 1)
    }else{
      j<-which(cumsum(probs)>(runif(1)*sum(probs)))[1]
    }
    gamma[j]<-1-gamma[j]
    ## update inclusion estimators
    if(t>burn_in){
      est_inclusion_probs<-est_inclusion_probs+gamma/T
    }
    ## store indices sequence for output analysis
    indices_sequence[,t]<-c(j_aux,j)
    gamma_1[t]<-gamma[1]
    gamma_2[t]<-gamma[2]
  }
  est_inclusion_probs_1<-cumsum(gamma_1[(burn_in+1):(burn_in+T)])/c(1:T)
  est_inclusion_probs_2<-cumsum(gamma_2[(burn_in+1):(burn_in+T)])/c(1:T)
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma,
              est_inclusion_probs_1=est_inclusion_probs_1,
              est_inclusion_probs_2=est_inclusion_probs_2))
}

#### AUXILIARY FUNCTIONS FOR SAMPLERS ######
###
single_full_cond<-function(j,gamma,hyper_par=NULL,stored=NULL){
  # hyper_par is a list containing y, X, prior_p_incl
  # and XtX, ytX, yty
  # y needs to be n x 1 matrix
  included<-which(gamma==1)
  p_gamma<-sum(gamma)
  n<-hyper_par$n
  p<-length(gamma)
  h<-hyper_par$prior_p_incl
  c<-hyper_par$c
  yty<-hyper_par$yty
  ytX<-hyper_par$ytX
  XtX<-hyper_par$XtX
  if(p_gamma==0){
    S_old<-yty
    S_new<-yty-c/(1+c)*ytX[j]^2/XtX[j,j]
    post_ratio<-h/(1-h)*(S_old/S_new)^(n/2)/(sqrt(1+c))
    fc_j<-1/(1+post_ratio)# fc is the full cond prob of the current value of gamma[j]
  }else{
    ytX_gamma<-matrix(ytX[included],nrow = 1,ncol = p_gamma)
    ytXFXty<-sum((solve(t(chol(XtX[included,included])))%*%t(ytX_gamma))^2)
    S_old<-yty-c/(1+c)*ytXFXty
    F<-solve(XtX[included,included])
    if(gamma[j]==0){
      xjtX<-matrix(XtX[j,included],nrow=1,ncol=length(included))
      d<-c(1/(XtX[j,j]-xjtX%*%F%*%t(xjtX)))
      ytXFXtxj<-ytX_gamma%*%F%*%t(xjtX)
      S_new<-yty-c/(1+c)*(ytXFXty+d*(ytXFXtxj-ytX[j])^2)
      post_ratio<-h/(1-h)*(S_old/S_new)^(n/2)/(sqrt(1+c))
      fc_j<-1/(1+post_ratio)# fc is the full cond prob of the current value of gamma[j]
    }
    if(gamma[j]==1){
      if(p_gamma==1){
        S_new<-yty
      }else{
        pos_rem<-which(included==j)
        included_new<-included[-pos_rem]
        ytX_gamma_new<-matrix(ytX[included_new],nrow = 1,ncol = p_gamma-1)
         F_new<-F[-pos_rem,-pos_rem]-(as.matrix(F[-pos_rem,pos_rem],nrow=p_gamma-1)%*%F[pos_rem,-pos_rem])/F[pos_rem,pos_rem]
        S_new<-yty-c/(1+c)*ytX_gamma_new%*%F_new%*%t(ytX_gamma_new)
      }
      post_ratio<-(1-h)/h*(S_old/S_new)^(n/2)*(sqrt(1+c))
      fc_j<-1/(1+post_ratio)# fc_j is the full cond prob of the current value of gamma[j]
    }
  }
  output<-list(fc_j=fc_j,stored=NULL)
  return(output)
}

full_cond<-function(gamma,hyper_par=NULL,stored=NULL){
    # hyper_par is a list containing y, X, prior_p_incl
  # and XtX, ytX, yty
  # y needs to be n x 1 matrix
  fc<-rep(NA,length(gamma))
  included<-which(gamma==1)
  p_gamma<-sum(gamma)
  n<-hyper_par$n
  p<-length(gamma)
  h<-hyper_par$prior_p_incl
  c<-hyper_par$c
  yty<-hyper_par$yty
  ytX<-hyper_par$ytX
  XtX<-hyper_par$XtX
  if(p_gamma==0){
    S_old<-yty
    for (j in 1:p){
      # computing S_new using inverse matrix or cholesky
      # X_gamma_new<-hyper_par$X[,j]
      # ytX_gamma_new<-ytX[j]
      S_new<-yty-c/(1+c)*ytX[j]^2/XtX[j,j]
      post_ratio<-h/(1-h)*(S_old/S_new)^(n/2)/(sqrt(1+c))
      fc[j]<-1/(1+post_ratio)# fc is the full cond prob of the current value of gamma[j]
    }
  }else{
    ytX_gamma<-matrix(ytX[included],nrow = 1,ncol = p_gamma)
    ytXFXty<-sum((solve(t(chol(XtX[included,included])))%*%t(ytX_gamma))^2)
    S_old<-yty-c/(1+c)*ytXFXty
    F<-solve(XtX[included,included])
    L<-chol(F)
    ## COMPUTE FIRST AS IF THEY WERE ALL gamma[j]=0
    d_vec<-1/(diag(XtX)-rowSums((XtX[,included]%*%t(L))^2))
    ytXFXtxj_vec<-matrix(ytX[included],nrow=1)%*%F%*%XtX[included,]
    S_new_vec<-yty-c/(1+c)*(ytXFXty+d_vec*(ytXFXtxj_vec-ytX)^2)
    post_ratio_vec<-h/(1-h)*(S_old/S_new_vec)^(n/2)/(sqrt(1+c))
    fc<-1/(1+post_ratio_vec)
    ## THEN MODIFY THE ONES FOR WHICH gamma[j]=1; COULD BE IMPROVED...
    for (j in included){
      if(p_gamma==1){
        S_new<-yty
      }else{
        pos_rem<-which(included==j)
        included_new<-included[-pos_rem]
        ytX_gamma_new<-matrix(ytX[included_new],nrow = 1,ncol = p_gamma-1)
        # using updating
        F_new<-F[-pos_rem,-pos_rem]-(as.matrix(F[-pos_rem,pos_rem],nrow=p_gamma-1)%*%F[pos_rem,-pos_rem])/F[pos_rem,pos_rem]
        S_new<-yty-c/(1+c)*ytX_gamma_new%*%F_new%*%t(ytX_gamma_new)
      }
      post_ratio<-(1-h)/h*(S_old/S_new)^(n/2)*(sqrt(1+c))
      fc[j]<-1/(1+post_ratio)# fc is the full cond prob of the current value of gamma[j]
    }
  }
  output<-list(fc=fc,stored=NULL)
  return(output)
}

compute_flip_rates<-function(output_fc,gamma){
  # TGS version
  flip_rates<-rep(NA,length(gamma))
  flip_rates[which(gamma==0)]<-(1-output_fc$fc[which(gamma==0)])/output_fc$fc[which(gamma==0)]
  flip_rates[which(gamma==1)]<-1
  return(flip_rates)
}

#### FUNCTION TO SIMULATE DATA ##########
###
simulate_data<-function(n,p,c,SNR=2,scenario=1,collinearity=100){
  require(MASS)
  Sigma<-diag(1,p)
  if(scenario==1){# variables 1 and 2 strongly correlated
    rho<-0.99
    Sigma[1,2]<-Sigma[2,1]<-rho
    X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
    beta0 <- matrix(c(1,rep(0,p-1)),ncol=1)
  }
  if(scenario==2){#correlated design of Wang et al (2011) and Huang et al (2016); with different variance
    rho<-0.9
    for (i in 1:3){for(j in 1:3){if(i!=j)Sigma[i,j]<-rho}}
    for (i in 4:6){for(j in 4:6){if(i!=j)Sigma[i,j]<-rho}}
    require(MASS)
    X <- mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma)
    beta0<-matrix(c(3,3,-2,3,3,-2,rep(0,p-6)),ncol=1)
  }
  if(scenario==3){# uncorrelated_design
    X <- matrix(rnorm(n*p),nrow=n,ncol=p)
    beta0 <- matrix(c(2,-3,2,2,-3,3,-2,3,-2,3,rep(0,p-10)),ncol=1)
  }
  sigma<-1
  beta<-SNR*sqrt(log(p)/n)*beta0
  y <- X %*% beta + rnorm(n,mean = 0,sd = sigma^2)
  X<-t(t(X)-colMeans(X))
  y<-y-mean(y)
  ### precompute matrices needed for samplers and create list to output ###
  XtX<-t(X)%*%X
  ytX<-t(y)%*%X
  yty<-sum(y^2)
  hyper_par<-list(
    n=n,y=y,X=X,XtX=XtX, ytX=ytX, yty=yty,
    prior_p_incl=prior_p_incl,
    c=c
  )
  return(hyper_par)
}
