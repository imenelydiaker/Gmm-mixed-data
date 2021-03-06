% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KmeansInit.R
\name{KmeansInit}
\alias{KmeansInit}
\title{Kmeans initialization for clustering with EM}
\usage{
KmeansInit(z, Xqt, Xcat, K, N, hj)
}
\arguments{
\item{z}{A list}

\item{Xqt}{A matrix}

\item{Xcat}{A matrix}

\item{K}{A list}

\item{N}{A number}

\item{hj}{A list}
}
\value{
Model parameters : theta = list(mu, mat_cov, pk, alpha)
}
\description{
Kmeans initialization for clustering with EM
}
