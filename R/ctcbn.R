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

ctcbn <- function(path)
{
  x = .Call("ctcbn", suppressWarnings(normalizePath(path)))
  rows = strsplit(x, "\n")[[1]]
  names = strsplit(rows[2], "\t")
  values = strsplit(rows[3], "\t")

  df <- data.frame(row.names = unlist(names))
  df$values = unlist(values)

  return(t(df))
}
