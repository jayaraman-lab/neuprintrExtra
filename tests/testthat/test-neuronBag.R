PFLNeurons <- getTypesTable("PFL1")
PFLFull1 <- neuronBag(PFLNeurons[1,],slctROI="LAL(R)")

test_that("neuronBag construction, concatenation and retyping works",{
  expect_is(PFLFull1,"neuronBag")
  PFLFull2 <- neuronBag(PFLNeurons[1,],by.roi=FALSE)
  expect_is(PFLFull2,"neuronBag")
  PFLFull3 <- neuronBag(PFLNeurons[1,])
  expect_is(PFLFull3,"neuronBag")
  
  expect_is(cxRetyping(PFLFull1),"neuronBag")
  
  PFLCombo <- c(PFLFull1,PFLFull3)
  expect_is(PFLCombo,"neuronBag")
})

test_that("Empty neuronBags are consistent", {
  PFLEmpty <- neuronBag(PFLNeurons[1,],slctROI="MB(R)")
  PFLEmptyPlus <- neuronBag(PFLNeurons[1,],slctROI="MB(R)",computeKnownRatio=TRUE)
  expect_is(PFLEmpty,"neuronBag")
  expect_is(PFLEmptyPlus,"neuronBag")
  expect_is(c(PFLEmptyPlus,PFLFull1),"neuronBag")
})
