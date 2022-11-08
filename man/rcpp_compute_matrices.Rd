\name{rcpp_compute_matrices}
\alias{rcpp_compute_matrices}
\title{Computes alignment matrix and directions matrix}
\usage{
rcpp_compute_matrices(Am, Dm, seq1, seq2, gap, miss, match)
}
\arguments{
  \item{Am}{Alignment matrix}
  
  \item{Dm}{Directions matrix}
  
  \item{seq1}{vector of strings containing the first sequence}
  
  \item{seq2}{vector of strings containing the second sequence}
  
  \item{gap}{An integer representing the score for a gap in the alignment}
  
  \item{miss}{An integer representing the score for a mismatch in the alignment}
  
  \item{match}{An integer representing the score for a match in the alignment}
}
\value{
Alignment an directions matrices filled
}
\description{
This should only be a function used by the main package function globalSequenceAlignment.
}
