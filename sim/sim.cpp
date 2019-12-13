#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sim_pop_cpp(const int& n_days,
                      const int& init_pop_size,
                      const NumericVector& pvac,
                      const NumericVector& pflu,
                      const double& ve,
                      const int& lag) {
  const int vec_size = n_days - 1;
  IntegerVector day(vec_size), unvac_uninf_sus(vec_size),
    vac_uninf_lag(lag), vac_uninf_sus(vec_size), vac_uninf_imm(vec_size),
    unvac_inf(vec_size), vac_inf(vec_size);
  unvac_uninf_sus[0] = init_pop_size;
  day[0] = 1;
  for (int i = 1; i < vec_size; i++) {
    // A = unvac_uninf_sus
    // B = vac_uninf_lag
    // C = vac_uninf_sus
    // D = vac_uninf_imm
    // E = unvac_inf
    // F = vac_inf
    int A_to_B = pvac[i - 1] * unvac_uninf_sus[i - 1];
    int A_to_E = pflu[i - 1] * unvac_uninf_sus[i - 1];
    int B_to_D = ve * vac_uninf_lag[lag - 1];
    int B_to_C = (1 - ve) * vac_uninf_lag[lag - 1];

    day[i] = i + 1;
  }
  DataFrame pop = DataFrame::create(
    Named("day") = day,
  );
  return pop;
}
