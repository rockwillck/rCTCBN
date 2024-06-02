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
ctcbn <- function(path,
                  bootstrap_mode = 0,
                  bootstrap_samples = 0,
                  random_seed = 0,
                  sampling_rate = 1.0,
                  epsilon = 2,
                  num_drawn_samples = 0,
                  num_em_runs = 1
)
{
  x = .Call("ctcbn", suppressWarnings(normalizePath(path)), as.integer(bootstrap_mode), as.integer(bootstrap_samples), as.integer(random_seed), as.double(sampling_rate), as.double(epsilon), as.integer(num_drawn_samples), as.integer(num_em_runs))
  rows = strsplit(x, "\n")[[1]]
  names = strsplit(rows[2], "\t")
  values = strsplit(rows[3], "\t")

  df <- data.frame(row.names = unlist(names))
  df$values = unlist(values)

  return(t(df))
}
