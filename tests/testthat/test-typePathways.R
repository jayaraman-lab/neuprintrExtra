FB6Neurons <- getTypesTable("FB6E")
FB6Bag <- neuronBag(FB6Neurons,by.roi=FALSE,selfRef=TRUE,renaming=cxRetyping)

test_that("Pathway functions work",{
  expect_is(pathDirect <- get_type2typePath(FB6Neurons,FB6Neurons,by.roi=FALSE,n_steps=1:2,renaming=cxRetyping),"data.frame")
  expect_is(pathDirect2 <- get_type2typePath(FB6Neurons,FB6Neurons,by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping),"data.frame")
  expect_is(pathFromBag <- tableChain2path(FB6Bag$outputs,FB6Bag$inputs,type.to=cxRetyping(FB6Neurons)),"data.frame")
  expect_equal(pathDirect,pathFromBag)
  
})

test_that("Contralateral pathway completion works",{
  FB1ANeurons <- getTypesTable("FB1A")
  FS1BNeurons <- getTypesTable("FS1B")
  expect_is(pathDirectBoth <- get_type2typePath(FS1BNeurons,FB1ANeurons,ROI = c("SNP(R)","FBl1"),n_steps=1:2,addContraPaths = T,renaming=cxRetyping),"data.frame")
  expect_is(pathDirectBoth2 <- get_type2typePath(FS1BNeurons,FB1ANeurons,ROI = c("SNP(R)","FBl1"),n_steps=1:2,addContraPaths = T,thresholdPerROI = 20,renaming=cxRetyping),"data.frame")
})