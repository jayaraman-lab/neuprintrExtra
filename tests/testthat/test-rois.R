test_that("ROI outline methods match", {
  no_mesh <- neuprintr::neuprint_ROI_mesh("NO(R)")
  noOutline <- roiOutline(no_mesh,roiName="NO(R)")
  expect_equal(noOutline,roiOutline("NO(R)"))
})

test_that("ROI innervation functions work",{
  expect_is(BULTable <- getTypesInRoiTable("BU(L)",renaming=lateralize_types),"neuronBag")
  expect_true(nrow(BULTable$outputs)>50)
  expect_is(BULTypes <- typesInROI(BULTable,"BU(L)"),"data.frame")
  expect_true("ExR1_L" %in% BULTypes$type)
})