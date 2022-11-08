
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_compute_matrices(IntegerMatrix Am, StringMatrix Dm,
                           StringVector seq1, StringVector seq2,
                           int gap, int miss, int match) {

    int nrow = Am.nrow(), ncol = Am.ncol();

    for (int i = 1; i < nrow; i++) {
      for (int j = 1; j < ncol; j++) {
        int vertical_score = Am(i-1, j) + gap;
        int horizontal_score = Am(i, j-1) + gap;
        int diagonal_score = 0;
        if (seq1[j-1] == seq2[i-1]) {
          diagonal_score = Am(i-1, j-1) + match;
        }
        else {
          diagonal_score = Am(i-1, j-1) + miss;
        }

        IntegerVector score = {vertical_score, horizontal_score, diagonal_score};

        int max_score = max(score);

        Am(i, j) = max_score;

        if (max_score == diagonal_score) {
          Dm(i, j) = "D";
        }
        else if (max_score == horizontal_score) {
          Dm(i, j) = "H";
        }
        else {
          Dm(i, j) = "V";
        }
      }
    }

    List results;
    results["alignment"] = Am;
    results["directions"] = Dm;

    return results;
}
