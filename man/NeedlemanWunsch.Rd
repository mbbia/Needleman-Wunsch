% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NeedlemanWunsch.R
\docType{package}
\name{NeedlemanWunsch}
\alias{NeedlemanWunsch}
\title{Perform pairwise global sequence alignment}
\description{
This package contains a function that perform pairwise global 
 sequence alignment using the Needleman-Wunsch algorithm
}
\details{
NeedlemanWunsch

This package contains a function that perform pairwise global sequence alignment
using the Needleman-Wunsch algorithm.

In the package is contained only the function globalSequenceAlignment.
 It performs the pairwise global sequence alignment of the two sequences and
 the parameters given by the user.
}
\examples{
\dontrun{
# Use vignette(NeedlemanWunsch) for more information
# create two DNAString objects from strings
seq1 <- DNAString('GAATC')
seq2 <- DNAString('CATACG')
  
# load DNAString objects from FASTA file
seq_from_fasta <- readDNAStringSet(pathToFile, format = "fasta")[[1]]
   
# use the main global sequence alignment function
## with default score parameters
globalSequenceAlignment(seq1, seq2)
   
## specifying score parameters
globalSequenceAlignment(seq1, seq2, -4, -5, 10)
   
# assign the results to a variable in order to access all the resulting variables
results <- globalSequenceAlignment(seq1, seq2)
   
# to access all the resulting variables
results$alignment   # a DNAStringSet containing the alignment 
results$score       # the final score of the alignment
results$score_mat   # the matrix filled with the scores
results$dir_mat     # the matrix filled with the directions for the traceback
}
}
\references{
Needleman-Wunsch algorithm:
 Needleman SB, Wunsch CD. A general method applicable to the search for similarities
 in the amino acid sequence of two proteins. Journal of molecular biology.
 1970 Mar 28;48(3):443-53
}
\author{
Marco Bianchi
 Mainteiner: Marco Bianchi
 E-Mail: <marco30.bianchi@mail.polimi.it> or <mb.bianchi.998@gmail.com>
}
