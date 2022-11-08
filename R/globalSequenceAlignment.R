#' Sequence alignment
#' 
#' This function performs pairwise global sequence alignment of two sequences.
#' It uses the parameters gap, miss and match in order to compute the alignment
#' matrix and the directions matrix and then prints the resulting alignment
#' between the two sequences
#' 
#' @param seq1 A sequence as Biostrings::DNAString object
#' @param seq2 A sequence as Biostrings::DNAString object
#' @param gap An integer representing the score for a gap in the alignment (default = -1)
#' @param miss An integer representing the score for a mismatch in the alignment (default = -1)
#' @param match An integer representing the score for a match in the alignment (default = +1)
#' @return A list containing the final alignment, the score of the alignment and the score and directions matrices.
globalSequenceAlignment <- function(seq1, seq2, gap = -1, miss = -1, match = 1) {
  
  # error handling
  stopifnot("seq1 must be a DNAString object"=class(seq1) == "DNAString")
  stopifnot("seq2 must be a DNAString object"=class(seq2) == "DNAString")
  stopifnot("seq1 length must be at least 1"=length(seq1) >= 1)
  stopifnot("seq2 length must be at least 1"=length(seq2) >= 1)
  stopifnot("gap score must be an integer"=all.equal(gap, as.integer(gap)))
  stopifnot("miss score must be an integer"=all.equal(miss, as.integer(miss)))
  stopifnot("match score must be an integer"=all.equal(match, as.integer(match)))
  
  # splitting the sequences in order to have vectors of strings
  seq1_split <- unlist(strsplit(toString(seq1), ""))
  seq2_split <- unlist(strsplit(toString(seq2), ""))
  
  # computing the length of the two sequences
  len1 <- length(seq1_split)
  len2 <- length(seq2_split)
  
  # creating the alignment score matrix
  alignment_matrix <- matrix(0, ncol = len1+1, nrow = len2+1)
  # filling first row and column of the alignment score matrix
  for (i in 2:nrow(alignment_matrix)) {
    alignment_matrix[i,1] <- (alignment_matrix[i,1]+(i-1))*(gap)
  }
  for (j in 2:ncol(alignment_matrix)) {
    alignment_matrix[1,j] <- (alignment_matrix[1,j]+(j-1))*(gap)
  }
  
  # creating the directions matrix
  directions_matrix <- matrix("0", ncol = len1+1, nrow = len2+1)
  # filling first row and column of the directions matrix
  directions_matrix[1,] <- rep("H")
  directions_matrix[,1] <- rep("V")
  directions_matrix[1,1] <- "0"
  
  # computing the resulting score and directions matrices with an Rcpp function
  resulting_matrices <- rcpp_compute_matrices(alignment_matrix, directions_matrix,
                                              seq1_split, seq2_split, gap, miss, match)
  
  # saving the output of the rcpp_compute_matrices function
  alignment_matrix <- resulting_matrices$alignment
  directions_matrix <- resulting_matrices$directions
  
  # rename the matrices axis with the sequences nucleotides
  rownames(directions_matrix) <- c("-", seq2_split)
  colnames(directions_matrix) <- c("-", seq1_split)
  
  # initialize the vector containing the alignments
  alignment1 <- c()
  alignment2 <- c()
  
  # traceback fase of the algorithm
  while (directions_matrix[i,j] != "0") {
    if (directions_matrix[i,j] == "H") {
      alignment1 <- append(alignment1, colnames(directions_matrix)[j], after = 0)
      alignment2 <- append(alignment2, "-", after = 0)
      j = j - 1
      
    } else if (directions_matrix[i,j] == "V") {
      alignment1 <- append(alignment1, "-", after = 0)
      alignment2 <- append(alignment2, rownames(directions_matrix)[i], after = 0)
      i = i - 1
      
    } else {
      alignment1 <- append(alignment1, colnames(directions_matrix)[j], after = 0)
      alignment2 <- append(alignment2, rownames(directions_matrix)[i], after = 0)
      i = i - 1
      j = j - 1
    }
  }
  
  # convert vector of strings into strings
  alignment1 <- base::paste(alignment1, collapse = "")
  alignment2 <- base::paste(alignment2, collapse = "")

  # saving the final alignment as DNAStringSet object
  final_alignment <- DNAStringSet(c(alignment1,alignment2))
  names(final_alignment) <- c("Sequence 1", "Sequence 2")
  # saving the score of the alignment
  alignment_score <- alignment_matrix[nrow(alignment_matrix), ncol(alignment_matrix)]
  
  # creating the list of results to return
  results <- list(alignment = final_alignment, score = alignment_score, 
                  score_mat = alignment_matrix, dir_mat = directions_matrix)
  
  return(results)
  
}
