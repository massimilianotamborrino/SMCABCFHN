#' @useDynLib guidedABCFHN
#' @importFrom Rcpp sourceCpp
NULL

#'@rdname abc_distance
#'@title abc_distance
#'@description Compute the distance given the summaries
#'@param summaries vector of summaries
#'@param summ_weights vector with weights of the summaries computed from a pilot study
#'@param weight weights for the model-based summaries
#'@param type_sum choose between classic, model based and mixed
#'@return ABC distance for the given summaries
#'@export
abc_distance<-function(type_sum,summaries,summ_weights,weight) {
  if(type_sum=='classic'|type_sum=='mixed') return(sqrt(sum(summaries/diag(summ_weights)* summaries)))
  else if (type_sum=='model-based') return(summaries[1,1]+weight*summaries[1,2])
  }

