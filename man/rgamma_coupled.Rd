% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma_couplings.R
\name{rgamma_coupled}
\alias{rgamma_coupled}
\title{Sample from maximally coupled Gamma variables}
\usage{
rgamma_coupled(alpha1, alpha2, beta1, beta2)
}
\arguments{
\item{alpha1}{First shape}

\item{alpha2}{Second shape}

\item{beta1}{First rate}

\item{beta2}{Second rate}
}
\value{
A list with entry 'xy' for the pair of values, and boolean 'identical' indicating whether the two values
are identical.
}
\description{
Draws a pair of variables, respectively Gamma(alpha1, beta1) and Gamma(alpha2, beta2)
where the parametrization is that beta is the rate, i.e. the log-pdf of Gamma(alpha,beta) evaluated at x is
\deqn{\alpha * log(\beta) - lgamma(\alpha) + (\alpha-1) * log(x) - \beta x}
where \eqn{lgamma} stands for the logarithm of the Gamma function.
}
