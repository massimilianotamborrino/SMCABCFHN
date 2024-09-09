#' @useDynLib SMCABCFHN
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname abc_summaries
#'@title abc_summaries
#'@description Summaries
#'@param x Matrix with observations;
#'@param extra_summaries List with any extra coefficients which may be needed in the summaries
#'@param whichsummarymodelbased specify whether to consider IAEspectrum, IAE density, Wass density, Wass spectrum
#'@return Invariant spectral density and invariant density
#'@export
abc_summaries<-function(x,extra_summaries,whichsummarymodelbased) {
if(extra_summaries[[1]]=='model-based'|extra_summaries[[1]]=='mixed') return(FHN_abc_summaries_X(x,extra_summaries[[2]],extra_summaries[[3]]))
  else if (extra_summaries[[1]]=='classic') return(abc_summaries_classic(x))
  else if (extra_summaries[[1]]=='mixedsum') return(matrix(c(abc_summaries_classic(x),FHN_abc_summaries(extra_summaries[[2]],x,extra_summaries[[3]],
                                             extra_summaries[[4]],extra_summaries[[5]],extra_summaries[[6]],extra_summaries[[7]],extra_summaries[[8]],whichsummarymodelbased)),ncol=1))
    else if (extra_summaries[[1]]=='IAEWass') {
      return(FHN_abc_summaries(extra_summaries[[2]],x,extra_summaries[[3]],extra_summaries[[4]],extra_summaries[[5]],extra_summaries[[6]],extra_summaries[[7]],extra_summaries[[8]],whichsummarymodelbased));
    }
}

  #'@rdname FHN_abc_summaries_X
  #'@title FHN_abc_summaries_X
  #'@description Summaries for the reference data X, namely spectral density and invariant density.
  #'@param X V component of the observed data
  #'@param span_val span value for the smoothing periodogram
  #'@param Lsupport double check this
#'@return List with step size for the invariant spectral density and invariant density, values of the spectral density and invariant density
  #'@export
FHN_abc_summaries_X<-function(X,span_val,Lsupport){
  specX<-spectrum(X,log="no",span=span_val,plot=FALSE)
  spx<-specX$freq#*(1/hsim)
  stepP<-diff(spx)[1]
  startSupp<--5
  endSupp<-5
  densX<-density(X,n=Lsupport,from=startSupp,to=endSupp)
  stepD<-diff(densX$x)[1]
  return(list(stepP,stepD,specX$spec,densX$y))}

#'@rdname FHN_abc_summaries
#'@title FHN_abc_summaries
#'@description FHN Parallel - 4 Parameters - Return the spectral density and the invariant density. V component is observed
#'@param X V component of the observed data
#'@param Y V component of the simulated data
#'@param span_val span value for the smoothing periodogram
#'@param Lsupport double check this
#'@param stepP step to use in the calculation of the Integrate Absolute Error for the spectral density
#'@param stepD step to use in the calculation of the Integrate Absolute Error for the invariant density
#'@param specX spectral density of X
#'@param invDensX invariant density of X
#'@param whichsummarymodelbased specify whether to consider IAEspectrum, IAE density, Wass density, Wass spectrum
#'@return See description
#'@export
FHN_abc_summaries<-function(X,Y,span_val,Lsupport,stepP,stepD,specX,invDensX,whichsummarymodelbased){
  specY<-spectrum(Y,log="no",span=span_val,plot=FALSE)$spec
  #-----------------------------------------------------------------------
  #Invariant Density: Manhattan
  #startSupp<-min(min(X),min(Y))
  # endSupp<-max(max(X),max(Y))
  #sum of distance (Area) between periodogram and invariant density
  startSupp<--5
  endSupp<-5
  invDensY<-density(Y,kernel="gaussian",from=startSupp,to=endSupp,n=Lsupport,adjust=1)$y
  {if(sum(whichsummarymodelbased==c(1,2))==2){
  summaries<-c(sum(abs(stepP*(specX-specY))),
               sum(abs(stepD*(invDensX-invDensY))))
  }
  else{
  wass<-wasserstein1d(X,Y)
  wass_spectrum<-wasserstein1d(specX,specY)
  wass_diff<-wasserstein1d(diff(X),diff(Y))
  summaries<-c(sum(abs(stepP*(specX-specY))),
               sum(abs(stepD*(invDensX-invDensY))),wass_spectrum,wass,wass_diff)[whichsummarymodelbased]}}
  d_vec<-matrix(summaries,ncol=1)
  #-----------------------------------------------------------------------
  return(d_vec)
}

#'@rdname abc_summaries_classic
#'@title abc_summaries_classic
#'@description Summaries
#'@param x Matrix with observations;
#'@return Summaries
#'@export
abc_summaries_classic<-function(x){
  d<-diff(x)
matrix(c(mean(x),var(x),skewness(x),kurtosis(x),acf(x,lag.max=5,plot=F)$acf[-1],mean(d),var(d),skewness(d),kurtosis(d),acf(d,lag.max=5,plot=F)$acf[-1]),ncol=1)
}


