#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame create_edge_list_cpp(IntegerVector nverts, IntegerVector simplices, IntegerVector times) {
  // First, calculate the total number of edges to pre-allocate memory
  long long total_edges = 0;
  for (int i = 0; i < nverts.size(); ++i) {
    long long k = nverts[i];
    if (k > 1) {
      total_edges += k * (k - 1) / 2;
    }
  }

  // Pre-allocate vectors for the final DataFrame
  IntegerVector author1(total_edges);
  IntegerVector author2(total_edges);
  IntegerVector year_col(total_edges);

  long long edge_idx = 0;
  long long simplex_idx = 0;

  // Main loop to process each publication
  for (int i = 0; i < nverts.size(); ++i) {
    int num_authors = nverts[i];
    if (num_authors > 1) {
      int current_year = times[i];
      // Generate pairs for the current publication
      for (int j = 0; j < num_authors; ++j) {
        for (int l = j + 1; l < num_authors; ++l) {
          author1[edge_idx] = simplices[simplex_idx + j];
          author2[edge_idx] = simplices[simplex_idx + l];
          year_col[edge_idx] = current_year;
          edge_idx++;
        }
      }
    }
    simplex_idx += num_authors;
  }

  return DataFrame::create(
    Named("author_1") = author1,
    Named("author_2") = author2,
    Named("year") = year_col
  );
}