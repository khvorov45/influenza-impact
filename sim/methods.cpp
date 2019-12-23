#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame method1_cpp(const IntegerVector& cases,
                      const IntegerVector& vaccinations,
                      const NumericVector& ve,
                      const int& nsam) {
  int nt = cases.size();
  IntegerVector pops(nt), popn(nt), casen(nt), averted_method1(nt);
  NumericVector pvac(nt), pflu(nt), vc_lag(nt);
  for (int i = 0; i < nt; i++) pvac[i] = vaccinations[i] / double(nsam);
  vc_lag[0] = pvac[0] / 2;
  pops[0] = nsam * (1 - vc_lag[0] * ve[0]);
  pflu[0] = double(cases[0]) / pops[0];
  popn[0] = nsam - cases[0];
  casen[0] = pflu[0] * popn[0];
  averted_method1[0] = casen[0] - cases[0];
  for (int i = 1; i < nt; i++) {
    vc_lag[i] = (pvac[i] + pvac[i - 1]) / 2;
    pops[i] = (pops[i - 1] - cases[i - 1]) * (1 - vc_lag[i] * ve[i]);
    pflu[i] = double(cases[i]) / pops[i];
    popn[i] = popn[i - 1] - casen[i - 1];
    casen[i] = pflu[i] * popn[i];
    averted_method1[i] = casen[i] - cases[i];
  }
  DataFrame method1 = DataFrame::create(
    _["cases"] = cases,
    _["vaccinations"] = vaccinations,
    _["ve"] = ve,
    _["vc_lag"] = vc_lag,
    _["pops"] = pops,
    _["pflu"] = pflu,
    _["popn"] = popn,
    _["casen"] = casen,
    _["averted_method1"] = averted_method1
  );
  return method1;
}

// [[Rcpp::export]]
DataFrame method2_cpp(const IntegerVector& cases,
                      const IntegerVector& vaccinations,
                      const NumericVector& ve,
                      const int& nsam) {
  int nt = cases.size();
  IntegerVector pop(nt), pops(nt), popn(nt), casen(nt), vef(nt),
    averted_method2(nt);
  NumericVector pvac(nt), pflu(nt);
  for (int i = 0; i < nt; i++) pvac[i] = vaccinations[i] / double(nsam);
  vef[0] = nsam * pvac[0] * ve[0];
  pop[0] = nsam - cases[0];
  pops[0] = nsam - cases[0] - vef[0];
  pflu[0] = cases[0] / double(nsam);
  popn[0] = nsam;
  casen[0] = pflu[0] * popn[0];
  averted_method2[0] = casen[0] - cases[0];
  for (int i = 1; i < nt; i++) {
    vef[i] = pop[i - 1] * pvac[i] * ve[i];
    pop[i] = pop[i - 1] - cases[i];
    pops[i] = pops[i - 1] - cases[i] - vef[i];
    pflu[i] = double(cases[i]) / pops[i - 1];
    popn[i] = popn[i - 1] - casen[i - 1];
    casen[i] = pflu[i] * popn[i];
    averted_method2[i] = casen[i] - cases[i];
  }
  DataFrame method2 = DataFrame::create(
    _["cases"] = cases,
    _["vaccinations"] = vaccinations,
    _["ve"] = ve,
    _["vef"] = vef,
    _["pop"] = pop,
    _["pops"] = pops,
    _["pflu"] = pflu,
    _["popn"] = popn,
    _["casen"] = casen,
    _["averted_method2"] = averted_method2
  );
  return method2;
}

// [[Rcpp::export]]
DataFrame method3_cpp(const IntegerVector& cases,
                      const IntegerVector& vaccinations,
                      const NumericVector& ve,
                      const int& nsam) {
  int nt = cases.size();
  IntegerVector b(nt), A(nt), C(nt), D(nt), E(nt), F(nt),
    popn(nt), casen(nt), averted_method3(nt);
  NumericVector pvac(nt), pflu(nt);
  pflu[0] = cases[0] / double(nsam);
  pvac[0] = vaccinations[0] / double(nsam);
  b[0] = nsam * pvac[0];
  A[0] = nsam * (1 - pflu[0]) - b[0];
  C[0] = b[0] * (1 - ve[0]);
  D[0] = b[0] * ve[0];
  E[0] = nsam * pflu[0];
  F[0] = 0;
  casen[0] = nsam * pflu[0];
  popn[0] = nsam - casen[0];
  averted_method3[0] = casen[0] - cases[0];
  for (int i = 1; i < nt; i++) {
    pflu[i] = double(cases[i]) / (A[i - 1] + C[i - 1]);
    pvac[i] = double(vaccinations[i]) / (A[i - 1] + E[i - 1]);
    b[i] = A[i - 1] * pvac[i];
    A[i] = A[i - 1] * (1 - pflu[i]) - b[i];
    C[i] = C[i - 1] * (1 - pflu[i]) + b[i] * (1 - ve[i]);
    D[i] = D[i - 1] + b[i] * ve[i];
    E[i] = E[i - 1] * (1 - pvac[i]) + A[i - 1] * pflu[i];
    F[i] = F[i - 1] + C[i - 1] * pflu[i] + E[i - 1] * pvac[i];
    casen[i] = pflu[i] * popn[i - 1];
    popn[i] = popn[i - 1] - casen[i];
    averted_method3[i] = casen[i] - cases[i];
  }
  DataFrame method3 = DataFrame::create(
    _["cases"] = cases,
    _["vaccinations"] = vaccinations,
    _["ve"] = ve,
    _["b"] = b,
    _["A"] = A,
    _["C"] = C,
    _["D"] = D,
    _["E"] = E,
    _["F"] = F,
    _["pvac"] = pvac,
    _["pflu"] = pflu,
    _["popn"] = popn,
    _["casen"] = casen,
    _["averted_method3"] = averted_method3
  );
  return method3;
}

// [[Rcpp::export]]
DataFrame method5_cpp(const IntegerVector& cases,
                      const IntegerVector& vaccinations,
                      const NumericVector& ve,
                      const int& nsam) {
  int nt = cases.size();
  IntegerVector casen(nt), averted_method5(nt);
  NumericVector pvac_cum(nt);
  pvac_cum[0] = vaccinations[0] / double(nsam);
  for (int i = 1; i < nt; i++) {
    pvac_cum[i] = pvac_cum[i - 1] + vaccinations[i] / double(nsam);
  }
  for (int i = 0; i < nt; i++) {
    casen[i] = cases[i] / (1 - ve[i] * pvac_cum[i]);
    averted_method5[i] = casen[i] - cases[i];
  }
  DataFrame method5 = DataFrame::create(
    _["cases"] = cases,
    _["vaccinations"] = vaccinations,
    _["pvac_cum"] = pvac_cum,
    _["casen"] = casen,
    _["averted_method5"] = averted_method5
  );
  return method5;
}
