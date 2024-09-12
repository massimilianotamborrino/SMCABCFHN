#' @useDynLib SMCABCFHN
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname abc_distance
#'@title abc_distance
#'@description Compute the distance given the summaries
#'@param summaries vector of summaries
#'@param summ_weights vector with weights of the summaries computed from a pilot study
#'@param weight weights for the structure-based summaries
#'@param type_sum choose between canonical and structure-based
#'@return ABC distance for the given summaries
#'@export
abc_distance<-function(type_sum,summaries,summ_weights,weight) {
  if(type_sum=='canonical') return(sqrt(sum(summaries/diag(summ_weights)* summaries)))
  else if (type_sum=='structure-based') return(summaries[1,1]+weight*summaries[1,2])
  }

