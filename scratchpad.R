library(ctcbn)

bc = dataset(poset = read_poset("~/Desktop/test/BC")$sets, numMutations = read_poset("~/Desktop/test/BC")$mutations, patternOrLambda = read_pattern("~/Desktop/test/BC"))
testingBC = ctcbn(bc, num_drawn_samples = 200)
