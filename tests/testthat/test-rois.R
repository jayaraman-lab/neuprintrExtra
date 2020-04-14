test_that("ROI outline methods match", {
  no_mesh <- neuprintr::neuprint_ROI_mesh("NO(R)")
  noOutline <- roiOutline(no_mesh,roiName="NO(R)")
  expect_equal(noOutline,roiOutline("NO(R)"))
})
