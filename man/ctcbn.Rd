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
\item{filestem}{string; Filestem of poset}
\item{bootstrap_mode}{logical;
FALSE: re-estimate poset (default), TRUE: re-estimate parameters for a fixed poset}
\item{bootstrap_samples}{integer; Number of bootstrap samples (requires epsilon > 0, num_drawn_samples = 0)}
\item{random_seed}{integer; Random number generator seed}
\item{sampling_rate}{float; Sampling rate \lambda_s = 1.0}
\item{epsilon}{float; If in [0,1], the fraction of violations allowed per edge; If negative, the interval [0,0.5] will be sampled equidistantly with N points; the resulting posets are written to the path <filestem>; If omitted, the poset in filestem.poset will be used}
\item{num_drawn_samples}{integer; If > 0, the number of samples to draw from the model; If zero (default), the model will be learned from data}
\item{num_em_runs}{integer; number of EM runs}
\value{
A dataframe
}
\description{
Runs ctcbn
}
\examples{
ctcbn("foo")
ctcbn("foo", epsilon=0.05)
ctcbn("foo", epsilon=-1)
}
