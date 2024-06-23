library(ctcbn)

bc = Spock$new(
  poset = read_poset("~/Desktop/test/BC")$sets,
  numMutations = read_poset("~/Desktop/test/BC")$mutations,
  patternOrLambda = read_pattern("~/Desktop/test/BC")
)
testingBC = ctcbn(bc, bootstrap_mode = TRUE, bootstrap_samples = 1000)

