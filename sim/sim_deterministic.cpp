#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sim_pop_cpp(const int& n_days,
                      const int& init_pop_size,
                      const IntegerVector& nvac,
                      const IntegerVector& nflu_novac,
                      const double& ve,
                      const int& lag) {
  List pop;

  // Day counter
  IntegerVector day(n_days);
  for (int i = 0; i < n_days; i++) {
    day[i] = i + 1;
  }
  pop["day"] = day;

  // A compartment (non-vaccinated, non-cases, susceptible)
  IntegerVector A(n_days);
  A[0] = init_pop_size;
  pop["A"] = A;
  IntegerVector A_to_B0(n_days);
  pop["A_to_B0"] = A_to_B0;
  IntegerVector A_to_E(n_days);
  pop["A_to_E"] = A_to_E;

  // B compartments (vaccinated, non-cases, lagging)
  for (int i = 0; i <= lag; i++) {
    std::string Bi_name = "B" + std::to_string(i);
    pop[Bi_name] = IntegerVector(n_days);
    std::string Bi_to_f_name = Bi_name + "_to_F";
    pop[Bi_to_f_name] = IntegerVector(n_days);
    if (i == lag) {
      std::string Bi_to_C_name = Bi_name + "_to_C";
      pop[Bi_to_C_name] = IntegerVector(n_days);
      std::string Bi_to_D_name = Bi_name + "_to_D";
      pop[Bi_to_D_name] = IntegerVector(n_days);
    }
  }

  // C compartment (vaccinated, non-cases, susceptible)
  IntegerVector C(n_days);
  pop["C"] = C;
  IntegerVector C_to_F(n_days);
  pop["C_to_F"] = C_to_F;

  // D compartment (vaccinated, non-cases, immune)
  IntegerVector D(n_days);
  pop["D"] = D;

  // E compartment (non-vaccinated cases)
  IntegerVector E(n_days);
  pop["E"] = E;
  IntegerVector E_to_F(n_days);
  pop["E_to_F"] = E_to_F;

  // F compartment (vaccinated cases)
  IntegerVector F(n_days);
  pop["F"] = F;

  // Parameter values
  pop["nvac"] = nvac;
  pop["nflu_novac"] = nflu_novac;
  pop["ve"] = ve;

  pop = DataFrame(pop);
  IntegerVector nvac_vec = pop["nvac"];
  IntegerVector nflu_novac_vec = pop["nflu_novac"];
  NumericVector ve_vec = pop["ve"];
  NumericVector pvac(n_days), pflu(n_days);
  pop["pvac"] = pvac;
  pop["pflu"] = pflu;
  IntegerVector uninf_novac(n_days), nflu(n_days);
  uninf_novac[0] = A[0];
  pop["uninf_novac"] = uninf_novac;
  pop["nflu"] = nflu;

  // Cycle days
  for (int i = 1; i < n_days; i++) {

    // Initialise compartments
    A[i] = A[i - 1];
    C[i] = C[i - 1];
    D[i] = D[i - 1];
    E[i] = E[i - 1];
    F[i] = F[i - 1];

    // Work out current parameter values
    pvac[i] = double(nvac_vec[i]) / (A[i - 1] + E[i - 1]);
    pflu[i] = double(nflu_novac_vec[i]) / uninf_novac[i - 1];
    uninf_novac[i] = uninf_novac[i - 1] - nflu_novac_vec[i];

    // A compartment (non-vaccinated, non-cases, susceptible)
    A_to_B0[i] = A[i - 1] * pvac[i - 1];
    A_to_E[i] = A[i - 1] * pflu[i - 1];
    A[i] -= A_to_B0[i] + A_to_E[i];
    nflu[i] += A_to_E[i];

    // B compartments (vaccinated, non-cases, lagging)
    for (int j = 0; j <= lag; j++) {
      std::string Bj_name = "B" + std::to_string(j);
      IntegerVector Bj = pop[Bj_name];
      std::string Bj_to_F_name = Bj_name + "_to_F";
      IntegerVector Bj_to_F = pop[Bj_to_F_name];
      Bj_to_F[i] = Bj[i - 1] * pflu[i - 1];
      F[i] += Bj_to_F[i];
      nflu[i] += Bj_to_F[i];
      if (j == 0) {
        Bj[i] = A_to_B0[i];
      } else {
        std::string Bjm1_name = "B" + std::to_string(j - 1);
        IntegerVector Bjm1 = pop[Bjm1_name];
        std::string Bjm1_to_F_name = Bjm1_name + "_to_F";
        IntegerVector Bjm1_to_F = pop[Bjm1_to_F_name];
        Bj[i] = Bjm1[i - 1] - Bjm1_to_F[i];
      }
      if (j == lag) {
        std::string Bj_to_D_name = Bj_name + "_to_D";
        IntegerVector Bj_to_D = pop[Bj_to_D_name];
        Bj_to_D[i] = (Bj[i - 1] - Bj_to_F[i])  * ve_vec[i - 1];
        std::string Bj_to_C_name = Bj_name + "_to_C";
        IntegerVector Bj_to_C = pop[Bj_to_C_name];
        Bj_to_C[i] = Bj[i - 1] - Bj_to_D[i] - Bj_to_F[i];
        C[i] += Bj_to_C[i];
        D[i] += Bj_to_D[i];
      }
    }

    // C compartment (vaccinated, non-cases, susceptible)
    C_to_F[i] = C[i - 1] * pflu[i - 1];
    C[i] -= C_to_F[i];
    nflu[i] += C_to_F[i];

    // E compartment (non-vaccinated cases)
    E_to_F[i] = E[i - 1] * pvac[i - 1];
    E[i] += A_to_E[i] - E_to_F[i];

    // F compartment (vaccinated cases)
    F[i] += C_to_F[i] + E_to_F[i];
  }

  return pop;
}
