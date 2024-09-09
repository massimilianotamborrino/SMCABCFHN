#' @useDynLib guidedABCFHN
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname abc_sequential_varying_threshold_number_sim
#'@title abc_sequential_varying_threshold_number_sim
#'@description Guided and non-guided SIS-ABC and SMC-ABC for time-varying thresholds computed as percentiles of ACCEPTED distances with a fixed budget in terms of number of simulations
#'@param data observed data
#'@param extra Any additional coefficients/vector of coefficients which may be need to simulate from the model
#'@param extra_summaries Any additional coefficients/vector needed to compute the summaries
#'@param parmask set of 0/1 to see which parameters have to be estimated
#'@param ABCthreshold Initial value for the threshold for the guided approach
#'@param number_sim Total number of simulations, after which we stop the algorithm
#'@param summ_weights vector with weights of the summaries computed from a pilot study
#'@param numparticles number of particles to keep
#'@param alpha percentile to compute the ABC threshold out of all distances
#'@param sampling specify the chosen samplers. The options are blocked, blockedopt, standard, olcm, full_cond, full_cond_opt
#'@param attempt number of iteration
#'@param nu degree of freedom of the t Copula (so type_copula =2)
#'@param type_copula 0 for the multivariate Gaussian, 1 for the Gaussian copula, 2 for the t-copula
#'@param type_margin 1 for triangular, 2 for t scaled, 3 for logistic, 5 for uniform, 6 for gaussian
#'@param nu_marg degree of freedom of the t location-scale marginal (so type_marginal =2)
#'@param folder specify the folder where the results are going to be saved in
#'@param we weights to scale the model-based summaries
#'@param whichsummarymodelbased specify whether to consider IAEspectrum, IAE density, Wass density, Wass spectrum
#'@param whichprior 1 for uniform, 2 for lognormal priors
#'@param subsamplingby every how many simulated points the observation should be taken. The default=1, i.e., no subsampling
#'@param whichiter Specify from when which iteration you want the new sampler to be used
#'@export
#'
abc_sequential_varying_threshold_number_sim<- function (data, extra, extra_summaries, parmask, ABCthreshold, number_sim, summ_weights, numparticles, alpha, sampling,
                                                        attempt, nu, type_copula, type_margin, nu_marg, folder, we,whichsummarymodelbased, whichprior = 1,subsamplingby=1,whichiter=5)
{
  #% THE OBSERVED SUMMARIES - Code for the FHN
  type_sum<-extra_summaries[[1]]
  span_val<-extra_summaries[[2]]
  Lsupport<-extra_summaries[[3]]
  extra_summaries_ref<-extra_summaries
  summobs<-abc_summaries(data,extra_summaries_ref);
  { if(type_sum=='model-based'){extra_summaries_ref<-list('IAEWass',data,span_val,Lsupport,summobs[[1]],summobs[[2]],summobs[[3]],summobs[[4]])
  summobs<-abc_summaries(data,extra_summaries_ref,whichsummarymodelbased)}
    else if(type_sum=='mixed'){extra_summaries_ref<-list('mixedsum',data,span_val,Lsupport,summobs[[1]],summobs[[2]],summobs[[3]],summobs[[4]])
    summobs<-abc_summaries(data,extra_summaries_ref,whichsummarymodelbased)}}

  #####
  nfreepar <- sum(parmask); # % the number of parameters to be inferred
  ABCdraws <- matrix(0,nrow=nfreepar,ncol=numparticles); #% will store accepted parameters (for 1 iteration only, not all)
  nsummaries <- length(summobs); #% the number of summary statistics
  simsumm_accepted <- matrix(0,nrow=nsummaries,ncol=numparticles);
  distance_accepted <- matrix(0,nrow=1,ncol=numparticles); #% this is only for the olcm proposal
  lengthmodel<-length(data)
  #% initialization: t is the iteration counter
  t <- 1;
  numproposals0<-0
  tic()   #% This will  count the seconds required to obtain the desired number of accepted particles for each iteration
  RES<-foreach(success = 1:numparticles,.combine='rbind',.packages='guidedABCFHN') %dopar% {
    distance<-ABCthreshold+1
    numproposals <- 0;
    simsumm_all <- c();
    while(distance >= ABCthreshold & numproposals<number_sim){
      numproposals <- numproposals +1;
      theta <-  problemprior(-1,1,whichprior);# % propose from the prior
      bigtheta <- theta# param_unmask(theta,parmask,parbase); % in this problem bigtheta <- theta. Don't bother
      simdata <- model(bigtheta,extra); #// simulate from the model
      simdata<- simdata[seq(1,length(simdata),by=subsamplingby)]
      simsumm <- abc_summaries(simdata,extra_summaries_ref,whichsummarymodelbased);# //% summaries from simdata
      xc <- (t(simsumm)-t(summobs)); #// compute
      distance <- abc_distance(type_sum,xc,summ_weights,we)
      simsumm_all <-  cbind(simsumm_all,simsumm)
    }
    if(numproposals>=number_sim) {list(numproposals);stop('a particle got stucked');}
    list(numproposals,distance,theta,simsumm,simsumm_all)
  }
  eval_time <- toc()
  weights <- matrix(1,ncol=numparticles);
  normweights <- weights/sum(weights); # % normalised weights
  ess <- 1/sum(normweights^2); # % the Effective Sample Size

  numproposals<-sum(unlist(RES[,1]))
  distance_accepted<-matrix(as.numeric(unlist(RES[,2])),nrow=1,ncol=numparticles)
  ABCdraws<-matrix(as.numeric(unlist(RES[,3])),nrow=nfreepar,ncol=numparticles)
  simsumm_accepted<-matrix(as.numeric(unlist(RES[,4])),nrow=nsummaries,ncol=numparticles)
  simsumm_all<-matrix(as.numeric(unlist(RES[,5])),nrow=nsummaries)

  totnumproposals<- numproposals
  rm(RES)

  write.table(eval_time,file=sprintf('%s/evaltime_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
  write.table(ess,file=sprintf('%s/ess_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
  write.table(ABCthreshold,file=sprintf('%s/ABCthreshold_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
  write.table(ABCdraws,file=sprintf('%s/ABCdraws_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
  write.table(numproposals,file=sprintf('%s/numproposals_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
  write.table(normweights,file=sprintf('%s/normweights_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)

  all_accepted_thetasimsum <- rbind(ABCdraws,simsumm_accepted);  #% all accepted parameters and summaries
  #% efficient vectorized computation of weighted covariance matrix for the
  #% full vector of parameters and summaries
  all_cov_thetasimsum <- t(all_accepted_thetasimsum) - matrix(rep(normweights %*% t(all_accepted_thetasimsum), numparticles),nrow=numparticles,byrow=T);  # % subtract weighted mean
  all_cov_thetasimsum <- t(all_cov_thetasimsum) %*% (all_cov_thetasimsum * matrix(rep(t(normweights), nfreepar+nsummaries),nrow=numparticles));
  all_cov_thetasimsum <- all_cov_thetasimsum / (1-sum(normweights^2)); #% the weighted covariance matrix
  all_cov_thetasimsum <- 0.5 * (all_cov_thetasimsum + t(all_cov_thetasimsum));   #% ensure symmetry
  dimthetasimsum<-dim(all_cov_thetasimsum)[1]

  ABCthreshold_temp <- quantile(distance_accepted,alpha/100);
  ABCthreshold_old <- ABCthreshold;
  ABCthreshold_new <- min(ABCthreshold_temp,ABCthreshold_old)#
  if (ABCthreshold_new == ABCthreshold_old)  ABCthreshold_new <- 0.95*ABCthreshold_new;

  if(sampling=='blocked'| sampling=='hybrid'| sampling=='mixed'| sampling=='mixedolcm'){
    mean_theta <- ABCdraws %*% t(normweights);  #% weighted mean across all accepted parameters
    mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
    cov_theta <- all_cov_thetasimsum[1:nfreepar,1:nfreepar]; # % parameters covariance matrix
    cov_simsum <- all_cov_thetasimsum[(nfreepar+1):dimthetasimsum,(nfreepar+1):dimthetasimsum]; #% summaries  covariance
    cov_thetasimsum <- all_cov_thetasimsum[1:nfreepar,(nfreepar+1):dimthetasimsum]; #% covariance of (theta,simsum)
    cov_simsumtheta <- t(cov_thetasimsum);
    proposal_mean <- mean_theta + cov_thetasimsum %*% (mldivide(cov_simsum,summobs-mean_simsum)); # % mean of the BLOCKED proposal
    proposal_mean <-  t(proposal_mean);  #% must be a row vector when passed to mvnrnd()
    proposal_cov <- cov_theta - cov_thetasimsum %*% mldivide(cov_simsum, cov_simsumtheta);  #% covariance of the BLOCKED proposal
    proposal_cov <- (proposal_cov + t(proposal_cov))/ 2; #% standard trick to ensure simmetry
    if(isposdef(proposal_cov)==0) proposal_cov<-nearPD(proposal_cov,base.matrix=TRUE)$mat
  }
  else if(sampling=="blockedopt"){
    mean_theta <- ABCdraws %*% t(normweights); # % weighted mean across all accepted parameters
    mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
    cov_simsum <- all_cov_thetasimsum[(nfreepar+1):dimthetasimsum,(nfreepar+1):dimthetasimsum];
    cov_thetasimsum <- all_cov_thetasimsum[1:nfreepar,(nfreepar+1):dimthetasimsum];
    cov_simsum <- (cov_simsum+t(cov_simsum))/2;
    if(isposdef(cov_simsum)==0) cov_simsum<-nearPD(cov_simsum,base.matrix=TRUE)$mat
    proposal_mean <- mean_theta + cov_thetasimsum %*% (mldivide (cov_simsum,summobs-mean_simsum));
    cov_blockedopt <-matrix(0,nrow=nfreepar,ncol=nfreepar);
    id_blockedopt <- which(distance_accepted < ABCthreshold_new);
    N0 <- length(id_blockedopt);  #% number of particles that satisfy the newest threshold
    if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
    { print('There are no particles satisfying the blockedopt criterion. We must leave..')
      return(1)}
    weights_blockedopt <- weights[id_blockedopt];
    normweights_blockedopt <- weights_blockedopt/sum(weights_blockedopt);
    for (jj in 1:N0) cov_blockedopt <- cov_blockedopt + normweights_blockedopt[jj]*(ABCdraws[,id_blockedopt[jj]]-proposal_mean)%*%t(ABCdraws[,id_blockedopt[jj]]-proposal_mean);
    cov_blockedopt <- (cov_blockedopt+t(cov_blockedopt))/2;
    if(isposdef(cov_blockedopt)==0) cov_blockedopt<-nearPD(cov_blockedopt,base.matrix=TRUE)$mat  # This IF statement is useful if the proposal_cov is not definite positive
    proposal_mean <- t(proposal_mean)# % must be a row vector when passed to mvnrnd()
  }
  else if(sampling=="full_cond"){
    mean_theta <- ABCdraws %*% t(normweights);  #% weighted mean across all accepted parameters
    mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
    proposal_mean <-proposal_var<- matrix(0,nrow=nfreepar,ncol=numparticles);
    for (ii in 1:numparticles){
      for (par in 1:nfreepar){
        #% notice the ~=par operator means "select all entries in the array except those corresponding to index par"
        proposal_mean[par,ii]<- mean_theta[par] + all_cov_thetasimsum[par,-par] %*% (mldivide(all_cov_thetasimsum[-par,-par],rbind(matrix(ABCdraws[-par,ii],ncol=1),summobs) - rbind(matrix(mean_theta[-par],ncol=1),mean_simsum)));
        proposal_var[par,ii]<- all_cov_thetasimsum[par,par] - all_cov_thetasimsum[par,-par] %*% (mldivide(all_cov_thetasimsum[-par,-par],all_cov_thetasimsum[-par,par]));
      }}
  }
  else if(sampling=="full_cond_opt"){
    mean_theta <- ABCdraws %*% t(normweights);  #% weighted mean across all accepted parameters
    mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
    proposal_mean <-var_opt_all<- matrix(0,nrow=nfreepar,ncol=numparticles);
    id_opt <- which(distance_accepted < ABCthreshold_new);
    N0 <- length(id_opt);  #% number of particles that satisfy the newest threshold
    if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
    { print('There are no particles satisfying the blockedopt criterion. We must leave..')
      return(2)}
    weights_opt <- weights[id_opt];
    normweights_opt <- weights_opt/sum(weights_opt);
    all_cov_thetasimsum <- (all_cov_thetasimsum+t(all_cov_thetasimsum))/2;
    if(isposdef(all_cov_thetasimsum)==0) all_cov_thetasimsum<-nearPD(all_cov_thetasimsum,base.matrix=TRUE)$mat  # This IF statement is useful if the proposal_cov is not definite positive
    for (ii in 1:numparticles)
      for (par in 1:nfreepar){
        proposal_mean[par,ii]<- mean_theta[par] + all_cov_thetasimsum[par,-par] %*% (mldivide(all_cov_thetasimsum[-par,-par],rbind(matrix(ABCdraws[-par,ii],ncol=1),summobs) - rbind(matrix(mean_theta[-par],ncol=1),mean_simsum)));
        for (jj in 1:N0) var_opt_all[par,ii]<- var_opt_all[par,ii] + normweights_opt[jj]*(ABCdraws[par,id_opt[jj]]-proposal_mean[par,ii])^2;
      }
  }
  else if(sampling=="standard"){
    C <- t(ABCdraws) - matrix(rep(normweights %*% t(ABCdraws), numparticles), nrow=numparticles,byrow=T);#   % subtract weighted mean                                                  % Remove mean (which is, also, weighted)
    C <- t(C) %*% (C * matrix(rep(t(normweights),nfreepar),nrow=numparticles));       # % Weighted Covariance Matrix
    C <- C / (1-sum(normweights^2)); #% the weighted covariance matrix
    C <- 0.5 * (C + t(C));   #% ensure symmetry
    Sigma <- 2*C;
  }
  else if (sampling=="olcm"){
    cov_olcm <- matrix(0,nrow=nfreepar,ncol=nfreepar);
    id_olcm <- which(distance_accepted < ABCthreshold_new);
    N0 <- length(id_olcm);  #% number of particles that satisfy the newest threshold
    if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
    { print('There are no particles satisfying the OLCM criterion. We must leave..')
      return(3)}
    weights_olcm <- weights[id_olcm];
    normweights_olcm <- weights_olcm/sum(weights_olcm);
    cov_olcm_all <- c();
    for (ii in 1:numparticles)
    {
      for (jj in 1:N0) cov_olcm <- cov_olcm + normweights_olcm[jj]*(ABCdraws[,id_olcm[jj]]-ABCdraws[,ii])%*%t(ABCdraws[,id_olcm[jj]]-ABCdraws[,ii]);
      cov_olcm <- (cov_olcm+t(cov_olcm))/2;
      cov_olcm_all <- cbind(cov_olcm_all,cov_olcm);
    }
  }

  #% refresh weighting matrix for ABC summaries
  tsimsumm_all<-t(simsumm_all)
  summ_weights<-diag(sapply(1:nsummaries,FUN=function(II) mean(abs(tsimsumm_all[,II]-mean(tsimsumm_all[,II])))^2))

  samplingoriginal<-0
  if(sampling=='hybrid') {sampling<- "blocked"
  samplingoriginal<-"hybrid"}

  if(sampling=='mixed') {sampling<- "blocked"
  samplingoriginal<-"mixed"}

  if(sampling=='mixedolcm') {sampling<- "blocked"
  samplingoriginal<-"mixedolcm"
  }

  while (totnumproposals < number_sim){ # I stop when the acceptance rate goes below 1.5% for the 2nd consecutive time
    t <- t+1
    weights_old <- weights;
    simsumm_accepted <- matrix(0,nrow=nsummaries,ncol=numparticles);
    if(t==whichiter& samplingoriginal=='mixed') type_margin<- 1
    if(t==whichiter& samplingoriginal=='mixedolcm') {sampling<-'olcm'
    ABCdraws_old <- ABCdraws;
    id_olcm <- which(distance_accepted < ABCthreshold_new);
    N0 <- length(id_olcm);  #% number of particles that satisfy the newest threshold
    if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
    { print('There are no particles satisfying the OLCM criterion. We must leave..')
      return(4)}
    weights_olcm <- weights_old[id_olcm];
    normweights_olcm <- weights_olcm/sum(weights_olcm);
    cov_olcm <- matrix(0,nrow=nfreepar,ncol=nfreepar);
    cov_olcm_all <- c();
    for (ii in 1:numparticles)
    {
      for (jj in 1:N0) cov_olcm <- cov_olcm + normweights_olcm[jj]*(ABCdraws[,id_olcm[jj]]-ABCdraws[,ii])%*%t(ABCdraws[,id_olcm[jj]]-ABCdraws[,ii]);
      cov_olcm <- (cov_olcm+t(cov_olcm))/2;
      cov_olcm_all <- cbind(cov_olcm_all,cov_olcm);
    }
    }

    if(sampling=='full_cond'){proposal_mean_old <- proposal_mean;
    proposal_var_old <- proposal_var;}
    else if(sampling=='standard'){ABCdraws_old <- ABCdraws;
    Sigma_old <- Sigma}
    else if(sampling=='olcm'){ABCdraws_old <- ABCdraws;
    #% now find the particles that would have been accepted (if any) under the NEW threshold
    id_olcm <- which(distance_accepted < ABCthreshold_new);
    N0 <- length(id_olcm);  #% number of particles that satisfy the newest threshold
    if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
    { print('There are no particles satisfying the OLCM criterion. We must leave..')
      return(4)}
 #   distance_accepted <-   matrix(0,nrow=1,ncol=numparticles); #% this is only for the olcm proposal
    weights_olcm <- weights_old[id_olcm];
    normweights_olcm <- weights_olcm/sum(weights_olcm);
    }
    else if(sampling=='full_cond_opt'){
      ABCdraws_old <- ABCdraws;
      proposal_mean_old <- proposal_mean;
      # now find the particles that would have been accepted (if any) under the NEW threshold
      id_opt <- which(distance_accepted < ABCthreshold_new);
      N0 <- length(id_opt);  #% number of particles that satisfy the newest threshold
      if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
      { print('There are no particles satisfying the full-cond-opt criterion. We must leave..')
        return(5)}
#      distance_accepted <- matrix(0,nrow=1,ncol=numparticles);
      weights_opt <- weights_old[id_opt];
      normweights_opt <- weights_opt/sum(weights_opt);
    }

    tic()

    RES2<-foreach(success = 1:numparticles,.combine='rbind',.packages=c('guidedABCFHN','pracma')) %dopar% {
      distance<- ABCthreshold_new+1
      numproposals <- 0;
      numproposals0 <- 0;
      simsumm_all<- c();
      theta<-rep(0,nfreepar)

      while(distance>=ABCthreshold_new & (numproposals+totnumproposals)<number_sim){
        numproposals <- numproposals +1;
        if(sampling=='blocked')  {sampletheta<-sample_theta(proposal_mean,proposal_cov,nu,type_copula,type_margin,nu_marg);
        theta<-sampletheta[[1]]
        }
        else if(sampling=='blockedopt') {sampletheta<-sample_theta(proposal_mean,cov_blockedopt,nu,type_copula,type_margin,nu_marg);
        theta<-sampletheta[[1]]
        }
        else if(sampling=='full_cond'){
          index <- sample(1:length(normweights),size=1,prob=normweights);
          for (par in 1:nfreepar) theta[par] <- proposal_mean_old[par,index] + sqrt(proposal_var_old[par,index]) * rnorm(1);
        }
        else if(sampling=='full_cond_opt'){
          index <- sample(1:length(normweights),size=1,prob=normweights);
          # id_opt contains all the indeces of the particles from the previous generation that satify the
          # current threshold. There are N0 of those.
          var_opt <- matrix(0,nrow=nfreepar,ncol=1);
          for (par in 1:nfreepar)
            for (jj in 1:N0) var_opt[par]<- var_opt[par] + normweights_opt[jj]*(ABCdraws_old[par,id_opt[jj]]-proposal_mean_old[par,index])^2;
          #          % the above variance is not "global" but is instead
          #         % specific for the sampled particle since it depends on proposal_mean_old(par,index)
          for (par in 1:nfreepar) theta[par]<- proposal_mean_old[par,index] + sqrt(var_opt[par]) * rnorm(1);
        }
        else if(sampling=='standard'){
          index <- sample(1:length(normweights),size=1,prob=normweights);
          theta <- rmvn(1,ABCdraws_old[,index],Sigma_old)
        }
        else if(sampling=='olcm'){
          index <- sample(1:length(normweights),size=1,prob=normweights);
          #          % id_olcm contains all the indeces of the particles from the previous generation that satisfy the
          #         % current threshold. There are N0 of those.
          cov_olcm <-matrix(0, nrow=nfreepar,ncol=nfreepar);
          for(jj in 1:N0) cov_olcm <- cov_olcm + normweights_olcm[jj]*(ABCdraws_old[,id_olcm[jj]]-ABCdraws_old[,index])%*%t(ABCdraws_old[,id_olcm[jj]]-ABCdraws_old[,index]);
          cov_olcm <- (cov_olcm+t(cov_olcm))/2;
          if(isposdef(cov_olcm)==0) cov_olcm<-nearPD(cov_olcm,base.matrix=TRUE)$mat  # This IF statement is useful if the proposal_cov is not definite positive
          #    % the above covariance is not "global" but is instead specific for the sampled particle
          theta <- rmvn(1,ABCdraws_old[,index],cov_olcm)
        }
        prior <- problemprior(theta,0,whichprior); #% evaluate prior
        if(prior==0) numproposals0<-numproposals0+1
        if(prior!=0){
          bigtheta <- theta#param_unmask(theta,parmask,parbase); % for this problem bigtheta <- theta
          simdata <- model(bigtheta,extra); #// simulate from the model
          simdata<- simdata[seq(1,length(simdata),by=subsamplingby)]
          simsumm <- abc_summaries(simdata,extra_summaries_ref,whichsummarymodelbased);# //% summaries from simdata
          xc <- (t(simsumm)-t(summobs)); #// compute
          distance <- abc_distance(type_sum,xc,summ_weights,we)
          simsumm_all <-  cbind(simsumm_all,simsumm);
        }#
      }
      if(numproposals>=number_sim) {list(numproposals);stop(print('a particle got stucked'));}
      # Now distance<ABCthreshold
      #  simsumm_accepted[,success]<-  simsumm; TO GET
      if(sampling=='blocked') {if (type_copula==0) {dens<- dmvn(theta,proposal_mean,proposal_cov);}
        else {dens <- dMvdc(theta,sampletheta[[2]])}}
      else if(sampling=='blockedopt') {
        if (type_copula==0) {dens<- dmvn(theta,proposal_mean,cov_blockedopt);}
        else {dens <- dMvdc(theta,sampletheta[[2]])}}
      else if(sampling=='full_cond'){
        dens <- 0;
        for(ii in 1:numparticles){
          denominator_prod <-1;
          for (par in 1:nfreepar) denominator_prod <- denominator_prod *dnorm(theta[par],proposal_mean_old[par,ii],sqrt(proposal_var_old[par,ii]));
          dens <- dens + weights_old[ii]*denominator_prod;}}
      else if(sampling=='full_cond_opt'){
        dens <- 0;
        for(ii in 1:numparticles){
          denominator_prod <-1;
          for (par in 1:nfreepar) denominator_prod <- denominator_prod *dnorm(theta[par],proposal_mean_old[par,ii],sqrt(var_opt_all[par,ii]));
          dens <- dens + weights_old[ii]*denominator_prod;}}
      else if(sampling=='standard'){
        dens <-0;
        for (ii in 1:numparticles) dens <- dens + weights_old[ii]*dmvn(theta,ABCdraws_old[,ii],Sigma_old);
      }
      else if(sampling=='olcm'){
        dens <-0;
        for (ii in 1:numparticles) dens <- dens + weights_old[ii]*dmvn(theta,ABCdraws_old[,ii],cov_olcm_all[,((ii-1)*nfreepar+1):(ii*nfreepar)]);
      }
      #  weights[success] <- prior /dens; TO GET
      #  ABCdraws[,success] <- theta; TO GET
      if(dens==0) return(print('error'))
      list(numproposals,distance,theta,simsumm,simsumm_all,prior/dens,numproposals0)
    }
    eval_time <- toc()

    numproposals<-sum(unlist(RES2[,1]))
    distance_accepted<-matrix(as.numeric(unlist(RES2[,2])),nrow=1,ncol=numparticles)
    ABCdraws<-matrix(as.numeric(unlist(RES2[,3])),nrow=nfreepar,ncol=numparticles)
    simsumm_accepted<-matrix(as.numeric(unlist(RES2[,4])),nrow=nsummaries,ncol=numparticles)
    simsumm_all<-matrix(as.numeric(unlist(RES2[,5])),nrow=nsummaries)
    weights<-matrix(as.numeric(unlist(RES2[,6])),ncol=numparticles);
    numproposals0<-sum(unlist(RES2[,7]))
    rm(RES2)

    totnumproposals<- totnumproposals+numproposals

    normweights <- weights/sum(weights); # % normalised weights
    ess <- 1/sum(normweights^2); # % the Effective Sample Size

    write.table(eval_time,file=sprintf('%s/evaltime_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
    write.table(ess,file=sprintf('%s/ess_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
    write.table(ABCthreshold_new,file=sprintf('%s/ABCthreshold_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
    write.table(ABCdraws,file=sprintf('%s/ABCdraws_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
    write.table(numproposals,file=sprintf('%s/numproposals_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
    write.table(numproposals0,file=sprintf('%s/numproposals0_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
    write.table(normweights,file=sprintf('%s/normweights_stage%d_attempt%d_copula%d_marg%d_nu%d_numarg%#.1f.txt',folder,t,attempt,type_copula,type_margin,nu,nu_marg),row.names = FALSE,col.names = FALSE)
    #   % BEFORE MOVING TO NEXT ITERATION, UPDATE THE COVARIANCES AND MEANS

    all_accepted_thetasimsum <- rbind(ABCdraws,simsumm_accepted);  #% all accepted parameters and summaries

    #% efficient vectorized computation of weighted covariance matrix for the
    #% full vector of parameters and summaries
    all_cov_thetasimsum <- t(all_accepted_thetasimsum) - matrix(rep(normweights %*% t(all_accepted_thetasimsum), numparticles),nrow=numparticles,byrow=T);  # % subtract weighted mean
    all_cov_thetasimsum <- t(all_cov_thetasimsum) %*% (all_cov_thetasimsum * matrix(rep(t(normweights), nfreepar+nsummaries),nrow=numparticles));
    all_cov_thetasimsum <- all_cov_thetasimsum / (1-sum(normweights^2)); #% the weighted covariance matrix
    all_cov_thetasimsum <- 0.5 * (all_cov_thetasimsum + t(all_cov_thetasimsum));   #% ensure symmetry

    ABCthreshold_temp <- quantile(distance_accepted,alpha/100);
    ABCthreshold_old <- ABCthreshold_new;
    ABCthreshold_new <- min(ABCthreshold_temp,ABCthreshold_old); #% ensure non-increasing thresholds
    if (ABCthreshold_new == ABCthreshold_old) {ABCthreshold_new <- 0.95*ABCthreshold_new;
    sprintf('Forced decrease of ABC threshold to: %f',ABCthreshold)
    }
    if(totnumproposals >= number_sim) return(ABCdraws)

    if(samplingoriginal=='hybrid') {sampling<- "blockedopt"
    samplingoriginal<-0}

    if(sampling=='blocked'){
      mean_theta <- ABCdraws %*% t(normweights);  #% weighted mean across all accepted parameters
      mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
      cov_theta <- all_cov_thetasimsum[1:nfreepar,1:nfreepar]; # % parameters covariance matrix
      cov_simsum <- all_cov_thetasimsum[(nfreepar+1):dimthetasimsum,(nfreepar+1):dimthetasimsum]; #% summaries  covariance
      cov_thetasimsum <- all_cov_thetasimsum[1:nfreepar,(nfreepar+1):dimthetasimsum]; #% covariance of (theta,simsum)
      cov_simsumtheta <- t(cov_thetasimsum);
      proposal_mean <- mean_theta + cov_thetasimsum %*% (mldivide(cov_simsum,summobs-mean_simsum)); # % mean of the BLOCKED proposal
      proposal_mean <-  t(proposal_mean);  #% must be a row vector when passed to mvnrnd()
      proposal_cov <- cov_theta - cov_thetasimsum %*% mldivide(cov_simsum, cov_simsumtheta);  #% covariance of the BLOCKED proposal
      proposal_cov <- (proposal_cov + t(proposal_cov))/ 2; #% standard trick to ensure simmetry
      if(isposdef(proposal_cov)==0) proposal_cov<-nearPD(proposal_cov,base.matrix=TRUE)$mat
    }
    else if(sampling=="blockedopt"){
      mean_theta <- ABCdraws %*% t(normweights); # % weighted mean across all accepted parameters
      mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
      cov_simsum <- all_cov_thetasimsum[(nfreepar+1):dimthetasimsum,(nfreepar+1):dimthetasimsum];
      cov_thetasimsum <- all_cov_thetasimsum[1:nfreepar,(nfreepar+1):dimthetasimsum];
      cov_simsum <- (cov_simsum+t(cov_simsum))/2;
      if(isposdef(cov_simsum)==0) cov_simsum<-nearPD(cov_simsum,base.matrix=TRUE)$mat
      proposal_mean <- mean_theta + cov_thetasimsum %*% mldivide(cov_simsum,summobs-mean_simsum);
      cov_blockedopt <- matrix(0,nrow=nfreepar,ncol=nfreepar);
      id_blockedopt<- which(distance_accepted < ABCthreshold_new);
      N0 <- length(id_blockedopt);  #% number of particles that satisfy the newest threshold
      if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
      {print('There are no particles satisfying the blockedopt criterion. We must leave...')
        return(5)}
      weights_blockedopt <- weights[id_blockedopt];
      normweights_blockedopt <- weights_blockedopt/sum(weights_blockedopt);
      for (jj in 1:N0) cov_blockedopt <- cov_blockedopt + normweights_blockedopt[jj]*(ABCdraws[,id_blockedopt[jj]]-proposal_mean)%*%t(ABCdraws[,id_blockedopt[jj]]-proposal_mean);
      cov_blockedopt <- (cov_blockedopt+t(cov_blockedopt))/2;
      if(isposdef(cov_blockedopt)==0) cov_blockedopt<-nearPD(cov_blockedopt,base.matrix=TRUE)$mat
      proposal_mean <- t(proposal_mean); }# % must be a row vector when passed to mvnrnd()

    else if(sampling=="full_cond"){
      mean_theta <- ABCdraws %*% t(normweights);  #% weighted mean across all accepted parameters
      mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
      proposal_mean <-proposal_var<- matrix(0,nrow=nfreepar,ncol=numparticles);
      all_cov_thetasimsum <- (all_cov_thetasimsum+t(all_cov_thetasimsum))/2;
      if(isposdef(all_cov_thetasimsum)==0) all_cov_thetasimsum<-nearPD(all_cov_thetasimsum,base.matrix=TRUE)$mat
      for (ii in 1:numparticles){
        for (par in 1:nfreepar){
          #% notice the ~=par operator means "select all entries in the array except those corresponding to index par"
          proposal_mean[par,ii]<- mean_theta[par] + all_cov_thetasimsum[par,-par] %*% (mldivide(all_cov_thetasimsum[-par,-par],rbind(matrix(ABCdraws[-par,ii],ncol=1),summobs) - rbind(matrix(mean_theta[-par],ncol=1),mean_simsum)));
          proposal_var[par,ii]<- all_cov_thetasimsum[par,par] - all_cov_thetasimsum[par,-par] %*% (mldivide(all_cov_thetasimsum[-par,-par],all_cov_thetasimsum[-par,par]));
        }}
    }
    else if(sampling=="full_cond_opt"){
      mean_theta <- ABCdraws %*% t(normweights);  #% weighted mean across all accepted parameters
      mean_simsum <- simsumm_accepted %*% t(normweights); #% weighted mean across all accepted summaries
      proposal_mean <-var_opt_all<- matrix(0,nrow=nfreepar,ncol=numparticles);
      id_opt <- which(distance_accepted < ABCthreshold_new);
      N0 <- length(id_opt);  #% number of particles that satisfy the newest threshold
      if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
      { print('There are no particles satisfying the blockedopt criterion. We must leave..')
        return(6)}
      weights_opt <- weights[id_opt];
      normweights_opt <- weights_opt/sum(weights_opt);
      all_cov_thetasimsum <- (all_cov_thetasimsum+t(all_cov_thetasimsum))/2;
      if(isposdef(all_cov_thetasimsum)==0) all_cov_thetasimsum<-nearPD(all_cov_thetasimsum,base.matrix=TRUE)$mat  # This IF statement is useful if the proposal_cov is not definite positive
      for (ii in 1:numparticles)
        for (par in 1:nfreepar){
          proposal_mean[par,ii]<- mean_theta[par] + all_cov_thetasimsum[par,-par] %*% (mldivide(all_cov_thetasimsum[-par,-par],rbind(matrix(ABCdraws[-par,ii],ncol=1),summobs) - rbind(matrix(mean_theta[-par],ncol=1),mean_simsum)));
          for (jj in 1:N0) var_opt_all[par,ii]<- var_opt_all[par,ii] + normweights_opt[jj]*(ABCdraws[par,id_opt[jj]]-proposal_mean[par,ii])^2;
        }
    }
    else if(sampling=="standard"){
      C <- t(ABCdraws) - matrix(rep(normweights %*% t(ABCdraws), numparticles), nrow=numparticles,byrow=T);#   % subtract weighted mean                                                  % Remove mean (which is, also, weighted)
      C <- t(C) %*% (C * matrix(rep(t(normweights),nfreepar),nrow=numparticles));       # % Weighted Covariance Matrix
      C <- C / (1-sum(normweights^2)); #% the weighted covariance matrix
      C <- 0.5 * (C + t(C));   #% ensure symmetry
      Sigma <- 2*C;
    }
    else if (sampling=="olcm"){
      cov_olcm <- matrix(0,nrow=nfreepar,ncol=nfreepar);
      id_olcm <- which(distance_accepted < ABCthreshold_new);
      N0 <- length(id_olcm);  #% number of particles that satisfy the newest threshold
      if (N0==0) # % in this sampling there are no particles satisfying the criterion above. We must exit
      { print('There are no particles satisfying the OLCM criterion. We must leave..')
        return(7)}
      weights_olcm <- weights[id_olcm];
      normweights_olcm <- weights_olcm/sum(weights_olcm);
      cov_olcm_all <- c();
      for (ii in 1:numparticles)
      {
        for (jj in 1:N0) cov_olcm <- cov_olcm + normweights_olcm[jj]*(ABCdraws[,id_olcm[jj]]-ABCdraws[,ii])%*%t(ABCdraws[,id_olcm[jj]]-ABCdraws[,ii]);
        cov_olcm <- (cov_olcm+t(cov_olcm))/2;
        cov_olcm_all <- cbind(cov_olcm_all,cov_olcm);
      }
    }

  }
  return(ABCdraws)
}
