# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

my_fun <- function(N, mu, Sigma) {
    .Call(`_SMCABCFHN_my_fun`, N, mu, Sigma)
}

FHN_prior3_ <- function(theta, draw) {
    .Call(`_SMCABCFHN_FHN_prior3_`, theta, draw)
}

FHN_prior2_ <- function(theta, draw) {
    .Call(`_SMCABCFHN_FHN_prior2_`, theta, draw)
}

FHN_prior_ <- function(theta, draw) {
    .Call(`_SMCABCFHN_FHN_prior_`, theta, draw)
}

FHN_model_ <- function(theta, delta, X0, N) {
    .Call(`_SMCABCFHN_FHN_model_`, theta, delta, X0, N)
}

