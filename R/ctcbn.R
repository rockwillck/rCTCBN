.onLoad <- function(lib, package)
{
  library.dynam("rCTCBN", package, lib)
}

ctcbn <- function(filestem,
                  bootstrap_mode = FALSE,
                  bootstrap_samples = 0,
                  random_seed = 0,
                  sampling_rate = 1.0,
                  epsilon = 2,
                  num_drawn_samples = 0,
                  num_em_runs = 1)
{
  x = .Call(
    "ctcbn",
    suppressWarnings(normalizePath(filestem)),
    as.integer(bootstrap_mode),
    as.integer(bootstrap_samples),
    as.integer(random_seed),
    as.double(sampling_rate),
    as.double(epsilon),
    as.integer(num_drawn_samples),
    as.integer(num_em_runs)
  )
  rows = strsplit(x, "\n")[[1]]
  names = strsplit(rows[2], "\t")
  values = strsplit(rows[3], "\t")

  df <- data.frame(row.names = unlist(names))
  df$values = unlist(values)

  return(as.data.frame(t(df)))
}
