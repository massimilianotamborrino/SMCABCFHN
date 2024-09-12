#' @useDynLib SMCABCFHN
#' @importFrom Rcpp sourceCpp evalCpp
NULL


#'@rdname model
#'@title model
#'@description Code for simulating using the splitting scheme from the FHN model
#'@param theta underlying parameter
#'@param extra Vector with time step delta, initial condition X0 and number of points N to simulate
#'@return Simulation of observations from all model, returned as matrix;
#'@export
model<-function(theta,extra) {
  FHN_model_(theta,extra[1],extra[2:3],extra[4])}

