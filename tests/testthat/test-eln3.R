test_that("eln3 is DNAString object", {
  seq <- NeedlemanWunsch::eln3
  expect_equal(class(seq)[[1]], "DNAString")
})