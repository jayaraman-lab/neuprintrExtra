PFLNeurons <- getTypesTable("PFL1")

test_that("Empty neuronBags are consistent", {
  PFLEmpty <- neuronBag(PFLNeurons[1,],slctROI="MB(R)")
  PFLEmptyPlus <- neuronBag(PFLNeurons[1,],slctROI="MB(R)",computeKnownRatio=TRUE)
  expect_is(PFLEmpty,"neuronBag")
  expect_is(PFLEmptyPlus,"neuronBag")
})
