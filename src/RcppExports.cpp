// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// disM
arma::mat disM(arma::mat xgrid, arma::vec ind);
RcppExport SEXP _BSPBSS_disM(SEXP xgridSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(disM(xgrid, ind));
    return rcpp_result_gen;
END_RCPP
}
// disM_full
arma::mat disM_full(arma::mat xgrid);
RcppExport SEXP _BSPBSS_disM_full(SEXP xgridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xgrid(xgridSEXP);
    rcpp_result_gen = Rcpp::wrap(disM_full(xgrid));
    return rcpp_result_gen;
END_RCPP
}
// samp_cov
arma::vec samp_cov(arma::mat X, arma::mat xgrid, arma::vec ind, int n);
RcppExport SEXP _BSPBSS_samp_cov(SEXP XSEXP, SEXP xgridSEXP, SEXP indSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(samp_cov(X, xgrid, ind, n));
    return rcpp_result_gen;
END_RCPP
}
// disM_full0
arma::mat disM_full0(arma::mat xgrid, double rho);
RcppExport SEXP _BSPBSS_disM_full0(SEXP xgridSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(disM_full0(xgrid, rho));
    return rcpp_result_gen;
END_RCPP
}
// cal_sumb
arma::mat cal_sumb(arma::mat& b, arma::mat& psi);
RcppExport SEXP _BSPBSS_cal_sumb(SEXP bSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_sumb(b, psi));
    return rcpp_result_gen;
END_RCPP
}
// cal_S
arma::mat cal_S(arma::mat& A, arma::mat& sumb, double zeta, double tau);
RcppExport SEXP _BSPBSS_cal_S(SEXP ASEXP, SEXP sumbSEXP, SEXP zetaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sumb(sumbSEXP);
    Rcpp::traits::input_parameter< double >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_S(A, sumb, zeta, tau));
    return rcpp_result_gen;
END_RCPP
}
// cal_S_new
arma::mat cal_S_new(arma::mat& A, arma::mat& sumb, double zeta);
RcppExport SEXP _BSPBSS_cal_S_new(SEXP ASEXP, SEXP sumbSEXP, SEXP zetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sumb(sumbSEXP);
    Rcpp::traits::input_parameter< double >::type zeta(zetaSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_S_new(A, sumb, zeta));
    return rcpp_result_gen;
END_RCPP
}
// cal_core
arma::mat cal_core(arma::mat& X, arma::mat& A, arma::mat& S);
RcppExport SEXP _BSPBSS_cal_core(SEXP XSEXP, SEXP ASEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_core(X, A, S));
    return rcpp_result_gen;
END_RCPP
}
// sum_core
double sum_core(arma::mat Xcore, arma::vec sigma);
RcppExport SEXP _BSPBSS_sum_core(SEXP XcoreSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xcore(XcoreSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_core(Xcore, sigma));
    return rcpp_result_gen;
END_RCPP
}
// dL_b_sub
arma::mat dL_b_sub(arma::mat& b, arma::mat& X, arma::mat& A, arma::vec& lambda, double& tau, arma::mat& psi, double& epsilon, double& zeta, arma::vec& sigma, int& sizep, int& sizen);
RcppExport SEXP _BSPBSS_dL_b_sub(SEXP bSEXP, SEXP XSEXP, SEXP ASEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP psiSEXP, SEXP epsilonSEXP, SEXP zetaSEXP, SEXP sigmaSEXP, SEXP sizepSEXP, SEXP sizenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int& >::type sizep(sizepSEXP);
    Rcpp::traits::input_parameter< int& >::type sizen(sizenSEXP);
    rcpp_result_gen = Rcpp::wrap(dL_b_sub(b, X, A, lambda, tau, psi, epsilon, zeta, sigma, sizep, sizen));
    return rcpp_result_gen;
END_RCPP
}
// GP_update_b_SGHMC
void GP_update_b_SGHMC(arma::mat& b, arma::mat& X, arma::mat& A, arma::mat& S, arma::vec& sigma, arma::vec& lambda, double& tau, arma::mat& psi, double& epsilon, arma::mat& sumb, double& zeta, arma::mat& X_core, double eta, double alpha, int sizep, int sizen, int m, int itr, arma::vec& nu);
RcppExport SEXP _BSPBSS_GP_update_b_SGHMC(SEXP bSEXP, SEXP XSEXP, SEXP ASEXP, SEXP SSEXP, SEXP sigmaSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP psiSEXP, SEXP epsilonSEXP, SEXP sumbSEXP, SEXP zetaSEXP, SEXP X_coreSEXP, SEXP etaSEXP, SEXP alphaSEXP, SEXP sizepSEXP, SEXP sizenSEXP, SEXP mSEXP, SEXP itrSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sumb(sumbSEXP);
    Rcpp::traits::input_parameter< double& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_core(X_coreSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type sizep(sizepSEXP);
    Rcpp::traits::input_parameter< int >::type sizen(sizenSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type itr(itrSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type nu(nuSEXP);
    GP_update_b_SGHMC(b, X, A, S, sigma, lambda, tau, psi, epsilon, sumb, zeta, X_core, eta, alpha, sizep, sizen, m, itr, nu);
    return R_NilValue;
END_RCPP
}
// GP_update_A
void GP_update_A(arma::mat& A, arma::vec& prior, arma::mat& X_core, arma::mat& X, arma::mat& S, arma::vec& sigma, int sizep);
RcppExport SEXP _BSPBSS_GP_update_A(SEXP ASEXP, SEXP priorSEXP, SEXP X_coreSEXP, SEXP XSEXP, SEXP SSEXP, SEXP sigmaSEXP, SEXP sizepSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_core(X_coreSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type sizep(sizepSEXP);
    GP_update_A(A, prior, X_core, X, S, sigma, sizep);
    return R_NilValue;
END_RCPP
}
// GP_update_sigma
void GP_update_sigma(arma::vec& sigma, arma::mat& X_core, arma::vec& prior_sigma);
RcppExport SEXP _BSPBSS_GP_update_sigma(SEXP sigmaSEXP, SEXP X_coreSEXP, SEXP prior_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_core(X_coreSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type prior_sigma(prior_sigmaSEXP);
    GP_update_sigma(sigma, X_core, prior_sigma);
    return R_NilValue;
END_RCPP
}
// log_p_zeta_Gaussian
double log_p_zeta_Gaussian(double zeta, arma::mat& X_core, arma::vec& sigma, arma::vec& prior);
RcppExport SEXP _BSPBSS_log_p_zeta_Gaussian(SEXP zetaSEXP, SEXP X_coreSEXP, SEXP sigmaSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_core(X_coreSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(log_p_zeta_Gaussian(zeta, X_core, sigma, prior));
    return rcpp_result_gen;
END_RCPP
}
// GP_update_zeta
void GP_update_zeta(double& zeta, double& tau, arma::mat& sumb, arma::mat& X_core, arma::vec& sigma, arma::mat& X, arma::mat& A, arma::mat& S, double stepsize, arma::vec prior);
RcppExport SEXP _BSPBSS_GP_update_zeta(SEXP zetaSEXP, SEXP tauSEXP, SEXP sumbSEXP, SEXP X_coreSEXP, SEXP sigmaSEXP, SEXP XSEXP, SEXP ASEXP, SEXP SSEXP, SEXP stepsizeSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sumb(sumbSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_core(X_coreSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    GP_update_zeta(zeta, tau, sumb, X_core, sigma, X, A, S, stepsize, prior);
    return R_NilValue;
END_RCPP
}
// GP_update_tau
void GP_update_tau(double& tau, double& zeta, arma::mat& sumb, arma::mat& X_core, arma::vec& sigma, arma::mat& X, arma::mat& A, arma::mat& S, double stepsize, arma::vec prior);
RcppExport SEXP _BSPBSS_GP_update_tau(SEXP tauSEXP, SEXP zetaSEXP, SEXP sumbSEXP, SEXP X_coreSEXP, SEXP sigmaSEXP, SEXP XSEXP, SEXP ASEXP, SEXP SSEXP, SEXP stepsizeSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double& >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sumb(sumbSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_core(X_coreSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior(priorSEXP);
    GP_update_tau(tau, zeta, sumb, X_core, sigma, X, A, S, stepsize, prior);
    return R_NilValue;
END_RCPP
}
// loglk
double loglk(arma::mat& X, arma::mat& A, arma::mat& S, arma::vec& sigma);
RcppExport SEXP _BSPBSS_loglk(SEXP XSEXP, SEXP ASEXP, SEXP SSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(loglk(X, A, S, sigma));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_bspbss_c
List mcmc_bspbss_c(arma::mat& X, arma::mat& A, arma::mat& b, double tau, arma::vec& sigma, double zeta, double subsample_n, double subsample_p, List prior, arma::mat& psi, arma::vec& lambda, double epsilon, double lr, double decay, int MClength, int burn_in, int thin, int show_step);
RcppExport SEXP _BSPBSS_mcmc_bspbss_c(SEXP XSEXP, SEXP ASEXP, SEXP bSEXP, SEXP tauSEXP, SEXP sigmaSEXP, SEXP zetaSEXP, SEXP subsample_nSEXP, SEXP subsample_pSEXP, SEXP priorSEXP, SEXP psiSEXP, SEXP lambdaSEXP, SEXP epsilonSEXP, SEXP lrSEXP, SEXP decaySEXP, SEXP MClengthSEXP, SEXP burn_inSEXP, SEXP thinSEXP, SEXP show_stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type zeta(zetaSEXP);
    Rcpp::traits::input_parameter< double >::type subsample_n(subsample_nSEXP);
    Rcpp::traits::input_parameter< double >::type subsample_p(subsample_pSEXP);
    Rcpp::traits::input_parameter< List >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type lr(lrSEXP);
    Rcpp::traits::input_parameter< double >::type decay(decaySEXP);
    Rcpp::traits::input_parameter< int >::type MClength(MClengthSEXP);
    Rcpp::traits::input_parameter< int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type show_step(show_stepSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_bspbss_c(X, A, b, tau, sigma, zeta, subsample_n, subsample_p, prior, psi, lambda, epsilon, lr, decay, MClength, burn_in, thin, show_step));
    return rcpp_result_gen;
END_RCPP
}
// smoos
NumericMatrix smoos(NumericMatrix S, IntegerMatrix xgrid, double smooth);
RcppExport SEXP _BSPBSS_smoos(SEXP SSEXP, SEXP xgridSEXP, SEXP smoothSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< double >::type smooth(smoothSEXP);
    rcpp_result_gen = Rcpp::wrap(smoos(S, xgrid, smooth));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BSPBSS_disM", (DL_FUNC) &_BSPBSS_disM, 2},
    {"_BSPBSS_disM_full", (DL_FUNC) &_BSPBSS_disM_full, 1},
    {"_BSPBSS_samp_cov", (DL_FUNC) &_BSPBSS_samp_cov, 4},
    {"_BSPBSS_disM_full0", (DL_FUNC) &_BSPBSS_disM_full0, 2},
    {"_BSPBSS_cal_sumb", (DL_FUNC) &_BSPBSS_cal_sumb, 2},
    {"_BSPBSS_cal_S", (DL_FUNC) &_BSPBSS_cal_S, 4},
    {"_BSPBSS_cal_S_new", (DL_FUNC) &_BSPBSS_cal_S_new, 3},
    {"_BSPBSS_cal_core", (DL_FUNC) &_BSPBSS_cal_core, 3},
    {"_BSPBSS_sum_core", (DL_FUNC) &_BSPBSS_sum_core, 2},
    {"_BSPBSS_dL_b_sub", (DL_FUNC) &_BSPBSS_dL_b_sub, 11},
    {"_BSPBSS_GP_update_b_SGHMC", (DL_FUNC) &_BSPBSS_GP_update_b_SGHMC, 19},
    {"_BSPBSS_GP_update_A", (DL_FUNC) &_BSPBSS_GP_update_A, 7},
    {"_BSPBSS_GP_update_sigma", (DL_FUNC) &_BSPBSS_GP_update_sigma, 3},
    {"_BSPBSS_log_p_zeta_Gaussian", (DL_FUNC) &_BSPBSS_log_p_zeta_Gaussian, 4},
    {"_BSPBSS_GP_update_zeta", (DL_FUNC) &_BSPBSS_GP_update_zeta, 10},
    {"_BSPBSS_GP_update_tau", (DL_FUNC) &_BSPBSS_GP_update_tau, 10},
    {"_BSPBSS_loglk", (DL_FUNC) &_BSPBSS_loglk, 4},
    {"_BSPBSS_mcmc_bspbss_c", (DL_FUNC) &_BSPBSS_mcmc_bspbss_c, 18},
    {"_BSPBSS_smoos", (DL_FUNC) &_BSPBSS_smoos, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_BSPBSS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
