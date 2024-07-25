library(rCTCBN)

bc = Spock$new(
  poset = matrix(),
  numMutations = 10,
  patternOrLambda = read_pattern("examples/BC")
)
testingBC = ctcbn(bc)


.Call("hcbn", "/Users/williamchoi-kim/Desktop/CTCBN_R/examples/BC")
