FB6Neurons <- getTypesTable("FB6E")
FB6Bag <- neuronBag(FB6Neurons,by.roi=FALSE,selfRef=TRUE,renaming=cxRetyping)

test_that("Pathway functions work",{
  expect_is(pathDirect <- get_type2typePath(FB6Neurons,FB6Neurons,by.roi=FALSE,n_steps=1:2,renaming=cxRetyping),"data.frame")
  expect_is(pathDirect2 <- get_type2typePath(FB6Neurons,FB6Neurons,by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping),"data.frame")
  expect_is(pathDirect3 <- get_type2typePath(FB6Neurons,FB6Neurons,by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping,ROI = "FB"),"data.frame")
  expect_is(pathDirect4 <- get_type2typePath(FB6Neurons,FB6Neurons,by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping,ROI = list("combo"=c("FB","SNP(R)"))),"data.frame")
  expect_is(pathFromBag <- tableChain2path(FB6Bag$outputs,FB6Bag$inputs,type.to=cxRetyping(FB6Neurons)),"data.frame")
  expect_equal(pathDirect,pathFromBag)
  
  openBag <- get_type2typePath_raw(FB6Neurons[1,],by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping,ROI = list("combo"=c("FB","SNP(R)")))
  expect_true(all(openBag[[2]]$type.from %in% openBag[[1]]$type.to))
  
  openBag2 <- get_type2typePath_raw(type.to=FB6Neurons[1,],by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping,ROI = list("combo"=c("FB","SNP(R)")))
  expect_true(all(openBag2[[1]]$type.to %in% openBag2[[2]]$type.from))
  
  openBag3 <- get_type2typePath_raw(FB6Neurons[1,],by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping,ROI = list("combo"=c("FB","SNP(R)")),thresholdPerROI = 20)
  expect_true(all(openBag3[[2]]$type.from %in% openBag3[[1]]$type.to))
  
  openBag4 <- get_type2typePath_raw(type.to=FB6Neurons[1,],by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping,ROI = list("combo"=c("FB","SNP(R)")),thresholdPerROI = 20)
  expect_true(all(openBag3[[1]]$type.to %in% openBag3[[2]]$type.from))
  
  openBag5 <- get_type2typePath_raw(type.to=FB6Neurons[2,],by.roi=FALSE,n_steps=1:2,stat=c("weightRelative","outputContribution"),renaming=cxRetyping,ROI = list("combo(R)"=c("SMP(R)")),thresholdPerROI = 20,addContraPaths = T)
  expect_true(all(openBag5[[1]]$type.to %in% openBag5[[2]]$type.from))
  
})

test_that("Contralateral pathway completion works",{
  FB1ANeurons <- getTypesTable("FB1A")
  FS1BNeurons <- getTypesTable("FS1B")
  expect_is(pathDirectBoth <- get_type2typePath(FS1BNeurons,FB1ANeurons,ROI = c("SNP(R)","FBl1"),n_steps=1:2,addContraPaths = T,renaming=cxRetyping),"data.frame")
  expect_is(pathDirectBoth2 <- get_type2typePath(FS1BNeurons,FB1ANeurons,ROI = c("SNP(R)","FBl1"),n_steps=1:2,addContraPaths = T,thresholdPerROI = 20,renaming=cxRetyping),"data.frame")
})