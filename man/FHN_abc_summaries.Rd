% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_summaries.r
\name{FHN_abc_summaries}
\alias{FHN_abc_summaries}
\title{FHN_abc_summaries}
\usage{
FHN_abc_summaries(
  X,
  Y,
  span_val,
  Lsupport,
  stepP,
  stepD,
  specX,
  invDensX,
  whichsummarymodelbased
)
}
\arguments{
\item{X}{V component of the observed data}

\item{Y}{V component of the simulated data}

\item{span_val}{span value for the smoothing periodogram}

\item{Lsupport}{double check this}

\item{stepP}{step to use in the calculation of the Integrate Absolute Error for the spectral density}

\item{stepD}{step to use in the calculation of the Integrate Absolute Error for the invariant density}

\item{specX}{spectral density of X}

\item{invDensX}{invariant density of X}

\item{whichsummarymodelbased}{specify whether to consider IAEspectrum, IAE density, Wass density, Wass spectrum}
}
\value{
See description
}
\description{
FHN Parallel - 4 Parameters - Return the spectral density and the invariant density. V component is observed
}
