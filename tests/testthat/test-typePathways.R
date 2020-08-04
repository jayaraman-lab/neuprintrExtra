FB6Neurons <- cxRetyping(getTypesTable("FB6E"))
FB6Bag <- cxRetyping(neuronBag(FB6Neurons[1,],by.roi=FALSE,selfRef=TRUE))

test_that("Pathway functions work",{
  expect_is(pathDirect <- get_type2typePath(FB6Neurons[1,],FB6Neurons[1,],by.roi=FALSE,n_steps=1:2),"data.frame")
  expect_is(pathFromBag <- tableChain2path(FB6Bag$outputs,FB6Bag$inputs,type.to=FB6Neurons[1,]),"data.frame")
  expect_equal(pathDirect,pathFromBag)
  
})

test_that("Contralateral pathway completion works",{
  FB1ANeurons <- cxRetyping(getTypesTable("FB1A"))
  FS1BNeurons <- cxRetyping(getTypesTable("FS1B"))
  expect_is(pathDirectBoth <- get_type2typePath(FS1BNeurons,FB1ANeurons,ROI = c("SNP(R)","CX"),n_steps=1:2,addContraPaths = T),"data.frame")
  expect_is(pathDirectBoth2 <- get_type2typePath(FS1BNeurons,FB1ANeurons,ROI = c("SNP(R)","CX"),n_steps=1:2,addContraPaths = T,thresholdPerROI = 20),"data.frame")
})