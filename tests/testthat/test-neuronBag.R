context("Neuron bags")

PFLNeurons <- getTypesTable("PFL1")
PFLFull1 <- neuronBag(PFLNeurons[1,],slctROI="LAL(R)")

test_that("neuronBag construction, concatenation and retyping works",{
  expect_is(PFLFull1,"neuronBag")
  PFLFull2 <- neuronBag(PFLNeurons[1,],by.roi=FALSE)
  expect_is(PFLFull2,"neuronBag")
  PFLFull3 <- neuronBag(PFLNeurons[1,])
  expect_is(PFLFull3,"neuronBag")
  PFLFull4 <- combineRois(PFLFull3,c("FB","LAL(R)"),"bigRegion")
  expect_is(PFLFull4,"neuronBag")
  
  expect_is(PFLFull1Ret <- cxRetyping(PFLFull1),"neuronBag")
  expect_is(PFLFull1RetAgain <- neuronBag(PFLNeurons[1,],slctROI="LAL(R)",renaming=cxRetyping),"neuronBag")
  expect_equal(nrow(PFLFull1RetAgain$inputs),nrow(PFLFull1Ret$inputs))
  expect_is(PFLRefd <- neuronBag(PFLNeurons[1,],computeKnownRatio = T),"neuronBag")
  expect_is(PFLRefdO <- neuronBag(PFLNeurons[1,],computeKnownRatio = T,omitInputs=T),"neuronBag")
  expect_is(PFLRefdI <- neuronBag(PFLNeurons[1,],computeKnownRatio = T,omitOutputs=T),"neuronBag")
  expect_true("ref" %in% names(PFLRefd))
  expect_equal(PFLRefd$outputs,PFLRefdO$outputs)
  expect_equal(PFLRefd$inputs,PFLRefdI$inputs)
  
  expect_true("ref" %in% names(PFLRefdO))
  expect_true("outputs_ref" %in% names(PFLRefdO$ref))
  expect_true("inputs_ref" %in% names(PFLRefdO$ref))
  
  PFLCombo <- c(PFLFull1,PFLFull3)
  expect_is(PFLCombo,"neuronBag")
  
  PFLComboRefd <- c(PFLRefd,PFLRefd)
  expect_is(PFLComboRefd,"neuronBag")
  expect_true("ref" %in% names(PFLComboRefd))
  
  PFLComboRefd <- combineRois(PFLComboRefd,c("LAL(R)","FB"),"combo")
  expect_is(PFLComboRefd,"neuronBag")
  expect_true("ref" %in% names(PFLComboRefd))
  
})

test_that("Empty neuronBags are consistent", {
  PFLEmpty <- neuronBag(PFLNeurons[1,],slctROI="MB(R)")
  PFLEmptyPlus <- neuronBag(PFLNeurons[1,],slctROI="MB(R)",computeKnownRatio=TRUE)
  expect_is(PFLEmpty,"neuronBag")
  expect_is(PFLEmptyPlus,"neuronBag")
  expect_is(c(PFLEmptyPlus,PFLFull1),"neuronBag")
})
