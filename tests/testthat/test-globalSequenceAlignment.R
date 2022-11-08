test_that("error DNAString seq1", {
  seq1 <- 'GAATC'
  seq2 <- DNAString("CATACG")
  expect_error(globalSequenceAlignment(seq1, seq2))
})

test_that("error DNAString seq2", {
  seq1 <- DNAString('GAATC')
  seq2 <- 'CATACG'
  expect_error(globalSequenceAlignment(seq1, seq2))
})

test_that("error DNAString length 1", {
  seq1 <- DNAString('')
  seq2 <- DNAString('CATACG')
  expect_error(globalSequenceAlignment(seq1, seq2))
})

test_that("error DNAString length 2", {
  seq1 <- DNAString('GAATC')
  seq2 <- DNAString('')
  expect_error(globalSequenceAlignment(seq1, seq2))
})

test_that("error gap variable", {
  seq1 <- DNAString('GAATC')
  seq2 <- DNAString('CATACG')
  gap <- -1.2
  miss <- -1
  match <- 1
  expect_error(globalSequenceAlignment(seq1, seq2, gap=gap, miss=miss, match=match))
})

test_that("error miss variable", {
  seq1 <- DNAString('GAATC')
  seq2 <- DNAString('CATACG')
  gap <- -1
  miss <- -1.4
  match <- 1
  expect_error(globalSequenceAlignment(seq1, seq2, gap=gap, miss=miss, match=match))
})

test_that("error match variable", {
  seq1 <- DNAString('GAATC')
  seq2 <- DNAString('CATACG')
  gap <- -1
  miss <- -1
  match <- 1.6
  expect_error(globalSequenceAlignment(seq1, seq2, gap=gap, miss=miss, match=match))
})

test_that("main function results 1", {
  seq1 <- DNAString("GAATC")
  seq2 <- DNAString("CATACG")
  observed <- globalSequenceAlignment(seq1,seq2)$alignment
  expected <- DNAStringSet(c('GA-ATC-', 'CATA-CG'))
  names(expected) <- c('Sequence 1', 'Sequence 2')
  expect_equal(observed, expected)
})

test_that("main function results 2", {
  seq1 <- DNAString("GAATC")
  seq2 <- DNAString("CATACG")
  observed <- globalSequenceAlignment(seq1,seq2)$score
  expected <- -1
  expect_equal(observed, expected)
})

test_that("main function results 3", {
  seq1 <- DNAString("GAATC")
  seq2 <- DNAString("CATACG")
  observed <- globalSequenceAlignment(seq1,seq2)$score_mat
  expected <- matrix(c(0,-1,-2,-3,-4,-5,
                       -1,-1,-2,-3,-4,-3,
                       -2,-2,0,-1,-2,-3,
                       -3,-3,-1,-1,0,-1,
                       -4,-4,-2,0,-1,-1,
                       -5,-5,-3,-1,-1,0,
                       -6,-4,-4,-2,-2,-1), nrow = 7, ncol = 6, byrow = TRUE)
  expect_equal(observed, expected)
})

test_that("main function results 3", {
  seq1 <- DNAString("GAATC")
  seq2 <- DNAString("CATACG")
  observed <- globalSequenceAlignment(seq1,seq2)$dir_mat
  expected <- matrix(c("0","H","H","H","H","H",
                       "V","D","D","D","D","D",
                       "V","D","D","D","H","H",
                       "V","D","V","D","D","H",
                       "V","D","D","D","H","D",
                       "V","D","V","V","D","D",
                       "V","D","V","V","D","V"), nrow = 7, ncol = 6, byrow = TRUE)
  rownames(expected) <- c("-",'C','A','T','A','C','G')
  colnames(expected) <- c("-",'G','A','A','T','C')
  expect_equal(observed, expected)
})
