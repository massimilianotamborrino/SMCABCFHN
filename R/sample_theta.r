#' @useDynLib guidedABCFHN
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname sample_theta
#'@title sample_theta
#'@param proposal_mean mean of the sampler
#'@param proposal_cov covariance matrix of the sampler
#'@param nu degree of freedom in case type_copula =2 (t student)
#'@param type_copula 0 for the multivariate Gaussian, 1 for the Gaussian copula, 2 for the t-copula, 3 for the Frank copula for 2-dim parameter vector
#'@param type_margin 1 for triangular, 2 for t scaled, 3 for logistic, 4 for lognormal, 5 for uniform, 6 for gaussian
#'@param nu_marg degree of freedom of the t location-scale marginal (so type_marginal =2)
#'@description Simulate 1 sample for theta
#'@return d-dimensional sampled value
#'@export
sample_theta<-function(proposal_mean,proposal_cov,nu,type_copula,type_margin,nu_marg){
#%%% Type_copula =0: Gaussian sampling
  #%%% Type copula = 1 Gaussian copula
  #%%% Type copula = 2 t copula
  #%%% type_margin = 1: triangular
  #%%% type_margin = 2  tlocation- scale
  #%%% type_margin = 3 logistic
  #%%% type_margin = 4 lognormal
  #%%% type_margin = 5 uniform
  #%%% type_margin = 6 gaussian

  margin<-c("triangle","lst","logis","lnorm","unif","norm")
  myMvd<-c();
  Dim<-dim(proposal_cov)[1];
  N<-1
  if (type_copula == 0)  {return(list(rmvn(N,proposal_mean,proposal_cov),myMvd));}
  else if (type_copula == 1) { rho<-cov2cor(proposal_cov);
  rho<- rho[lower.tri(rho)]
  myCop <-  normalCopula(rho,dim=Dim,dispstr="un")
  }
  else  if (type_copula== 2) {rho<-cov2cor(proposal_cov);
  rho<- rho[lower.tri(rho)]
  myCop <- tCopula(rho,df=nu,dim=Dim,dispstr="un")};
#  if (type_copula== 3){
#    rho<-cov2cor(proposal_cov)[1,2];
#  myCop<-frankCopula(BiCopTau2Par(5,rho, check.taus = TRUE),2);}
  # myMvd <- mvdc(copula=myCop, margins=c("gamma", "beta", "t"),
  #                 paramMargins=list(list(shape=2, scale=1),
  #                                   list(shape1=2, shape2=2),
  #                                   list(df=5)) )
    MARG<-list()
    if (type_margin ==  1){
      for (i in 1:Dim)  MARG[[i]]<-list(a=proposal_mean[i]-sqrt(6*proposal_cov[i,i]),b=proposal_mean[i]+sqrt(6*proposal_cov[i,i]),c=proposal_mean[i])
      myMvd<-mvdc(copula=myCop, margins=rep(margin[type_margin],Dim), paramMargins=MARG)
      theta <- rMvdc(N,myMvd)
      return(list(theta,myMvd))
      }
    if (type_margin == 2){
   for (i in 1:Dim) MARG[[i]]<-list(mu=proposal_mean[i],sigma=sqrt((nu_marg-2)*proposal_cov[i,i]/nu_marg),df=nu_marg);
   myMvd<- mvdc(copula=myCop, margins=rep(margin[type_margin],Dim), paramMargins=MARG)
   theta <- rMvdc(N,myMvd)
   return(list(theta,myMvd))}
    if (type_margin ==  3){
      for (i in 1:Dim) MARG[[i]]<-list(location=proposal_mean[i],scale=sqrt(3*proposal_cov[i,i])/pi);
      myMvd<-mvdc(copula=myCop, margins=rep(margin[type_margin],Dim), paramMargins=MARG)
        theta <- rMvdc(N,myMvd)
        return(list(theta,myMvd))}
   if (type_margin == 5 ){
      for (i in 1:Dim) MARG[[i]]<- list(min=proposal_mean[i]-sqrt(3*proposal_cov[i,i]),max=proposal_mean[i]+sqrt(3*proposal_cov[i,i]));
      myMvd<-mvdc(copula=myCop, margins=rep(margin[type_margin],Dim), paramMargins=MARG)
    theta <- rMvdc(N,myMvd)
    return(list(theta,myMvd))}
  if (type_margin == 4 ){
    for (i in 1:Dim) MARG[[i]]<- list(meanlog=log(proposal_mean[i]),sdlog=sqrt(log(0.5+0.5*sqrt(1+4*proposal_cov[i,i]/(proposal_mean[i]^2)))));
    myMvd<-mvdc(copula=myCop, margins=rep(margin[type_margin],Dim), paramMargins=MARG)
    theta <- rMvdc(N,myMvd)
    return(list(theta,myMvd))}

 if (type_margin == 6){
      for (i in 1:Dim) MARG[[i]]<- list(mean=proposal_mean[i],sd=sqrt(proposal_cov[i,i]))
      myMvd<-mvdc(copula=myCop, margins=rep(margin[type_margin],Dim),paramMargins=MARG)
        theta <- rMvdc(N,myMvd)
        return(list(theta,myMvd))}
    }
