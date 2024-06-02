\name{ctcbn}
\alias{ctcbn}
\title{Perform the ctcbn analysis.}
\usage{
ctcbn(path,
      bootstrap_mode = 0,
      bootstrap_samples = 0,
      random_seed = 0,
      sampling_rate = 1.0,
      epsilon = 2,
      num_drawn_samples = 0,
      num_em_runs = 1,
      x = 0,
      p = 0
      )
}
\arguments{
\item{path}{string; filestem of path}
\item{bootstrap_mode}{logical; determines whether or not bootstrapping is used}
\item{bootstrap_samples}{integer; number of bootstrap samples}
\item{random_seed}{integer; random seed}
\item{sampling_rate}{float; sampling rate (0.0 < n <  1.0)}
\item{epsilon}{float; convergence limit (0.0 < epsilon < 1.0)}
\item{num_drawn_samples}{integer; number of samples drawn}
\item{num_em_runs}{integer; number of EM runs}
\value{
A dataframe
}
\description{
Runs ctcbn
}
\examples{
ctcbn("CRC")
}
