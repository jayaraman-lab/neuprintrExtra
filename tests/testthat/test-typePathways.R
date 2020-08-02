FB6Neurons <- cxRetyping(getTypesTable("FB6E"))
FB6Bag <- cxRetyping(neuronBag(FB6Neurons[1,],by.roi=FALSE,selfRef=TRUE))

test_that("Pathway functions work",{
  expect_is(pathDirect <- get_type2typePath(FB6Neurons[1,],FB6Neurons[1,],by.roi=FALSE,n_steps=1:2),"data.frame")
  expect_is(pathFromBag <- tableChain2path(FB6Bag$outputs,FB6Bag$inputs,type.to=FB6Neurons[1,]),"data.frame")
  expect_equal(pathDirect,pathFromBag)
})