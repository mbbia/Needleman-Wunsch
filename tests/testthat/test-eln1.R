test_that("eln1 is DNAString object", {
  seq <- NeedlemanWunsch::eln1
  expect_equal(class(seq)[[1]], "DNAString")
})