#' @useDynLib SMCABCFHN
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname SMCABC_numbersim
#'@title SMCABC_numbersim
#'@description SMC-ABC for time-varying thresholds computed as percentiles of ACCEPTED distances with a fixed budget (number of simulations)
#'@param data observed data
#'@param extra Vector with time step delta, initial condition X0 and number of points N to simulate
#'@param extra_summaries Any additional coefficients/vector needed to compute the summaries
#'@param ABCthreshold Initial value for the threshold for the guided approach
#'@param number_sim Total number of simulations, after which we stop the algorithm
#'@param summ_weights vector with weights of the summaries computed from a pilot study
#'@param numparticles number of particles to keep
#'@param alpha percentile to compute the ABC threshold out of all distances
#'@param sampling specify the chosen samplers. The options are standard or olcm
#'@param attempt number of iteration
#'@param folder specify the folder where the results are going to be saved in
#'@param we weights to scale the structure-based summaries
#'@param whichsummarymodelbased specify whether to consider IAEspectrum, IAE density, Wass density, Wass spectrum
#'@param whichprior choose between 'unif','lognormal' and and 'exp'
#'@param subsamplingby every how many simulated points the observation should be taken. The default=1, i.e., no subsampling
#'@export
#'
SMCABC_numbersim<- function (data, extra, extra_summaries, ABCthreshold, number_sim, summ_weights, numparticles, alpha, sampling,
                                                        attempt, folder, we,whichsummarymodelbased, whichprior = 'unif',subsamplingby=1)
{
  #% THE OBSERVED SUMMARIES - Code for the FHN
  type_sum<-extra_summaries[[1]]
  span_val<-extra_summaries[[2]]
  Lsupport<-extra_summaries[[3]]
  extra_summaries_ref<-extra_summaries
  summobs<-abc_summaries(data,extra_summaries_ref);
  { if(type_sum=='structure-based'){extra_summaries_ref<-list('IAEWass',data,span_val,Lsupport,summobs[[1]],summobs[[2]],summobs[[3]],summobs[[4]])
  summobs<-abc_summaries(data,extra_summaries_ref,whichsummarymodelbased)}
    }

  #####
  nfreepar <- 4 # % the number of parameters to be inferred
  ABCdraws <- matrix(0,nrow=nfreepar,ncol=numparticles); #% will store accepted parameters (for 1 iteration only, not all)
  nsummaries <- length(summobs); #% the number of summary statistics
  distance_accepted <- matrix(0,nrow=1,ncol=numparticles); #% this is only for the olcm proposal
  lengthmodel<-length(data)
  #% initialization: t is the iteration counter
  t <- 1;
  numproposals0<-0
  numproposalsneg<-0
  tic()   #% This will  count the seconds required to obtain the desired number of accepted particles for each iteration
  RES<-foreach(success = 1:numparticles,.combine='rbind',.packages='SMCABCFHN') %dopar% {
    distance<-ABCthreshold+1
    numproposals <- 0;
    simsumm_all <- c();
    while(distance >= ABCthreshold & numproposals<number_sim){
      numproposals <- numproposals +1;
      theta <-  problemprior(-1,1,whichprior);# % propose from the prior
      bigtheta <- theta
      simdata <- model(bigtheta,extra); #// simulate from the model
      simdata<- simdata[seq(1,length(simdata),by=subsamplingby)]
      simsumm <- abc_summaries(simdata,extra_summaries_ref,whichsummarymodelbased);# //% summaries from simdata
      xc <- (t(simsumm)-t(summobs)); #// compute
      distance <- abc_distance(type_sum,xc,summ_weights,we)
      simsumm_all <-  cbind(simsumm_all,simsumm)
    }
    if(numproposals>=number_sim) {list(numproposals);stop('a particle got stucked');}
    list(numproposals,distance,theta,simsumm_all)
  }
  eval_time <- toc()
  weights <- matrix(1,ncol=numparticles);
  normweights <- weights/sum(weights); # % normalised weights
  ess <- 1/sum(normweights^2); # % the Effective Sample Size

  numproposals<-sum(unlist(RES[,1]))
  distance_accepted<-matrix(as.numeric(unlist(RES[,2])),nrow=1,ncol=numparticles)
  ABCdraws<-matrix(as.numeric(unlist(RES[,3])),nrow=nfreepar,ncol=numparticles)
  simsumm_all<-matrix(as.numeric(unlist(RES[,4])),nrow=nsummaries)

  #% refresh weighting matrix for ABC summaries
  tsimsumm_all<-t(simsumm_all)
  summ_weights<-diag(sapply(1:nsummaries,FUN=function(II) mean(abs(tsimsumm_all[,II]-mean(tsimsumm_all[,II])))^2))

  totnumproposals<- numproposals
  rm(RES)

  write.table(eval_time,file=sprintf('%s/evaltime_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(ess,file=sprintf('%s/ess_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(ABCthreshold,file=sprintf('%s/ABCthreshold_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(ABCdraws,file=sprintf('%s/ABCdraws_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(numproposals,file=sprintf('%s/numproposals_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
  write.table(normweights,file=sprintf('%s/normweights_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)

  ABCthreshold_temp <- quantile(distance_accepted,alpha/100);
  ABCthreshold_old <- ABCthreshold;
  ABCthreshold_new <- min(ABCthreshold_temp,ABCthreshold_old)#
  if (ABCthreshold_new == ABCthreshold_old)  ABCthreshold_new <- 0.95*ABCthreshold_new;

 if(sampling=="standard"){
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

  samplingoriginal<-0

  while (totnumproposals < number_sim){ # I stop when the acceptance rate goes below 1.5% for the 2nd consecutive time
    t <- t+1
    weights_old <- weights;
    if(sampling=='standard'){ABCdraws_old <- ABCdraws;
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

    tic()

    RES2<-foreach(success = 1:numparticles,.combine='rbind',.packages=c('SMCABCFHN','pracma')) %dopar% {
      distance<- ABCthreshold_new+1
      numproposals <- 0;
      numproposals0 <- 0;
      simsumm_all<- c();
      theta<-rep(0,nfreepar)

      while(distance>=ABCthreshold_new & (numproposals+totnumproposals)<number_sim){
        numproposals <- numproposals +1;
        if(sampling=='standard'){
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
        if(min(theta)<0) numproposalsneg<-numproposalsneg+1
        prior <- problemprior(theta,0,whichprior); #% evaluate prior
        if(prior==0) numproposals0<-numproposals0+1
        if(prior!=0){
          bigtheta <- theta
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
      if(sampling=='standard'){
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
      list(numproposals,distance,theta,simsumm_all,prior/dens,numproposals0,numproposalsneg)
    }
    eval_time <- toc()

    numproposals<-sum(unlist(RES2[,1]))
    distance_accepted<-matrix(as.numeric(unlist(RES2[,2])),nrow=1,ncol=numparticles)
    ABCdraws<-matrix(as.numeric(unlist(RES2[,3])),nrow=nfreepar,ncol=numparticles)
    simsumm_all<-matrix(as.numeric(unlist(RES2[,4])),nrow=nsummaries)
    weights<-matrix(as.numeric(unlist(RES2[,5])),ncol=numparticles);
    numproposals0<-sum(unlist(RES2[,6]))
    numproposalsneg<-sum(unlist(RES2[,7]))
    rm(RES2)

    totnumproposals<- totnumproposals+numproposals

    normweights <- weights/sum(weights); # % normalised weights
    ess <- 1/sum(normweights^2); # % the Effective Sample Size

    write.table(eval_time,file=sprintf('%s/evaltime_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(ess,file=sprintf('%s/ess_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(ABCthreshold_new,file=sprintf('%s/ABCthreshold_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(ABCdraws,file=sprintf('%s/ABCdraws_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(numproposals,file=sprintf('%s/numproposals_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(numproposals0,file=sprintf('%s/numproposals0_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(numproposalsneg,file=sprintf('%s/numproposalsneg_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)
    write.table(normweights,file=sprintf('%s/normweights_stage%d_attempt%d.txt',folder,t,attempt),row.names = FALSE,col.names = FALSE)

    ABCthreshold_temp <- quantile(distance_accepted,alpha/100);
    ABCthreshold_old <- ABCthreshold_new;
    ABCthreshold_new <- min(ABCthreshold_temp,ABCthreshold_old); #% ensure non-increasing thresholds
    if (ABCthreshold_new == ABCthreshold_old) {ABCthreshold_new <- 0.95*ABCthreshold_new;
    sprintf('Forced decrease of ABC threshold to: %f',ABCthreshold)
    }
    if(totnumproposals >= number_sim) return(ABCdraws)

    if(sampling=="standard"){
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
