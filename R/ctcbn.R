library(R6)

.onLoad <- function(lib, package)
{
  library.dynam("rCTCBN", package, lib)
}

read_pattern <- function(filestem) {
  lines <- readLines(paste(suppressWarnings(normalizePath(filestem)),
                           ".pat",
                           sep = ""))

  # Parse dimensions from the first line
  dimensions <- as.numeric(strsplit(lines[1], "\\s+")[[1]][2:3])
  num_rows <- dimensions[1]
  num_cols <- dimensions[2]

  # Read matrix data from subsequent lines
  matrix_data <- sapply(lines[2:length(lines)], function(line)
    as.numeric(strsplit(line, "\\s+")[[1]]))

  # Convert matrix data into a matrix
  matrix_data <- matrix(matrix_data, nrow = num_rows, byrow = FALSE)
  return(t(matrix_data))
}
read_poset <- function(filestem) {
  lines <- readLines(paste(suppressWarnings(normalizePath(filestem)), ".poset", sep =
                             ""))
  allSets = c()
  for (i in 3:length(lines) - 1) {
    allSets = c(allSets, as.numeric(unlist(
      strsplit(lines[[i]], " ")
    )[[1]]))
    allSets = c(allSets, as.numeric(unlist(
      strsplit(lines[[i]], " ")
    )[[2]]))
  }

  return(list(mutations = as.numeric(lines[[1]]), sets = matrix(allSets, ncol=2, byrow = TRUE)))
}
read_lambda <- function(filestem) {
  lines <- unlist(as.numeric(readLines(
    paste(suppressWarnings(normalizePath(filestem)), ".lambda", sep = "")
  )))
  return(matrix(lines, ncol = 1))
}

matrix_to_string <- function(mat) {
  # Get the dimensions of the matrix
  nrows <- nrow(mat)
  ncols <- ncol(mat)

  # Initialize an empty string to store the result
  result <- ""

  # Loop through each row
  for (i in 1:nrows) {
    res <- paste(mat[i, ], collapse = " ")
    # Loop through each column
    # for (j in 1:ncols) {
    #   # Append the matrix element to the result string
    #   result <- paste(result, mat[i, j], sep = " ")
    # }
    # Add a newline character at the end of each row (except the last row)
    if (i == 1) {
      result <- res

    } else if (i <= nrows) {
      result <- paste(result, res, sep = "\n")
    }
  }

  return(result)
}

temp_file = function(ext, fileC) {
  tf = tempfile("tempo", fileext = ext)[[1]]
  tfObj = file(tf)
  write(fileC, tfObj)
  close(tfObj)
  return(substr(tf, 0, nchar(tf) - nchar(ext)))
}

Spock = R6Class("Spock", list(
  poset = list(),
  numMutations = 0,
  patternOrLambda = matrix(),
  initialize = function (poset, numMutations, patternOrLambda) {
    stopifnot(is.matrix(poset))
    stopifnot(is.numeric(numMutations))
    stopifnot(is.matrix(patternOrLambda))
    self$poset = poset
    self$numMutations = numMutations
    self$patternOrLambda = patternOrLambda
  },
  getPoset = function() {
    output = ""
    if (nrow(self$poset) == 2) {
      for (i in 1:ncol(self$poset)) {
        output = paste(output, paste(self$poset[i, 1], self$poset[i, 2]), sep = "")
        output = paste(output, "\n")
      }
    }
    fileC = (paste(paste(
      as.character(self$numMutations), output, sep = "\n"
    ), "0", sep = ""))
    return(temp_file(ext = ".poset", fileC = fileC))
  },
  getSecond = function(n) {
    if (n == 0) {
      return(self$getPattern())
    } else {
      return(self$getLambda())
    }
  },
  getPattern = function() {
    fileC = (paste(
      paste(as.character(nrow(self$patternOrLambda)), as.character(ncol(self$patternOrLambda)), sep = " "),
      matrix_to_string(self$patternOrLambda),
      sep = "\n"
    ))

    return(temp_file(ext = ".pat", fileC = fileC))
  },
  getLambda = function() {
    output = ""
    for (i in 1:length(self$patternOrLambda)) {
      output = paste(output, self$patternOrLambda[[i]], sep = "")
      output = paste(output, "\n")
    }
    fileC = (output)
    return(temp_file(ext = ".lambda", fileC = fileC))
  }
))

filter_strings_by_start <- function(strings, start_substring) {
  filtered_strings <- grep(paste0("^", start_substring), strings, value = TRUE)
  return(filtered_strings)
}

ctcbn <- function(datasetObj,
                  bootstrap_mode = FALSE,
                  bootstrap_samples = 0,
                  random_seed = 0,
                  sampling_rate = 1.0,
                  epsilon = 2,
                  num_drawn_samples = 0,
                  num_em_runs = 1)
{
  outputStem = tempfile("output")
  x = .Call(
    "ctcbn",
    outputStem,
    datasetObj$getPoset(),
    datasetObj$getSecond(num_drawn_samples),
    as.integer(bootstrap_mode),
    as.integer(bootstrap_samples),
    as.integer(random_seed),
    as.double(sampling_rate),
    as.double(epsilon),
    as.integer(num_drawn_samples),
    as.integer(num_em_runs)
  )

  splitted = unlist(strsplit(x, "/"))

  outDir = paste(splitted[1:(length(splitted) - 1)], collapse = "/")
  outFiles = (filter_strings_by_start(list.files(outDir), splitted[[length(splitted)]]))

  outputList = list()
  for (file in outFiles) {
    file = paste(c(outDir, file), collapse = "/")
    if (endsWith(file, ".poset")) {
      outputList$poset = read_poset(substring(file, 1, nchar(file) - 6))
    }
    if (endsWith(file, ".pat")) {
      outputList$pattern = read_pattern(substring(file, 1, nchar(file) - 4))
    }
    if (endsWith(file, ".lambda")) {
      outputList$lambda = read_lambda(substring(file, 1, nchar(file) - 7))
    }
  }
  return(outputList)
}

hcbn <- function(datasetObj,
                  anneal = FALSE,
                  temp = 0,
                  annealing_steps = 0)
{
  outputStem = tempfile("output")
  x = .Call(
    "hcbn",
    outputStem,
    datasetObj$getPoset(),
    datasetObj$getSecond(num_drawn_samples),
    as.integer(anneal),
    as.double(temp),
    as.integer(annealing_steps)
  )
  
  splitted = unlist(strsplit(x, "/"))
  
  outDir = paste(splitted[1:(length(splitted) - 1)], collapse = "/")
  outFiles = (filter_strings_by_start(list.files(outDir), splitted[[length(splitted)]]))
  
  outputList = list()
  for (file in outFiles) {
    file = paste(c(outDir, file), collapse = "/")
    if (endsWith(file, ".poset")) {
      outputList$poset = read_poset(substring(file, 1, nchar(file) - 6))
    }
    if (endsWith(file, ".pat")) {
      outputList$pattern = read_pattern(substring(file, 1, nchar(file) - 4))
    }
    if (endsWith(file, ".lambda")) {
      outputList$lambda = read_lambda(substring(file, 1, nchar(file) - 7))
    }
  }
  return(outputList)
}
