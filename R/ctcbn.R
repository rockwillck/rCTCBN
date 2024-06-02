# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
.onLoad <- function(lib, package)
{
  library.dynam("rCTCBN", package, lib )
}

ctcbn <- function(path, bootstrap_mode, bootstrap_samples, random_seed, sampling_rate, epsilon=2, num_drawn_samples, num_em_runs, x, p)
{
  x = .Call("ctcbn", suppressWarnings(normalizePath(path)), as.integer(bootstrap_mode), as.integer(bootstrap_samples), as.integer(random_seed), as.double(sampling_rate), as.double(epsilon), as.integer(num_drawn_samples), as.integer(num_em_runs), as.logical(x), as.logical(p))
  rows = strsplit(x, "\n")[[1]]
  names = strsplit(rows[2], "\t")
  values = strsplit(rows[3], "\t")

  df <- data.frame(row.names = unlist(names))
  df$values = unlist(values)

  return(t(df))
}
