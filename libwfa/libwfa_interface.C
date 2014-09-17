#include "libwfa_interface.h"
#include "wf_analysis.h"

using namespace libwfa;

void analyse_opdm(const char *name, const char *desc,
    double *ddm_a, double *dm0_a, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;

    ab_matrix ddm(true), dm0(true);
    ddm.alpha() = arma::mat(ddm_a, nao, nao, false);
    dm0.alpha() = arma::mat(dm0_a, nao, nao, false);
    wf_analysis_static::analysis->analyse_opdm(std::cout,
            std::string(name), std::string(desc), ddm, dm0);
}


void analyse_opdm(const char *name, const char *desc,
    double *ddm_a, double *ddm_b, double *dm0_a, double *dm0_b, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;

    ab_matrix ddm(false), dm0(false);
    ddm.alpha() = arma::mat(ddm_a, nao, nao, false);
    ddm.beta()  = arma::mat(ddm_b, nao, nao, false);
    dm0.alpha() = arma::mat(dm0_a, nao, nao, false);
    dm0.beta()  = arma::mat(dm0_b, nao, nao, false);
    wf_analysis_static::analysis->analyse_opdm(std::cout,
            std::string(name), std::string(desc), ddm, dm0);
}


void analyse_opsdm(const char *name, const char *desc,
    double *dm_a, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;

    ab_matrix dm(true);
    dm.alpha() = arma::mat(dm_a, nao, nao, false);
    wf_analysis_static::analysis->analyse_opdm(std::cout,
            std::string(name), std::string(desc), dm);
}


void analyse_opsdm(const char *name, const char *desc,
    double *dm_a, double *dm_b, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;

    ab_matrix dm(false);
    dm.alpha() = arma::mat(dm_a, nao, nao, false);
    dm.beta()  = arma::mat(dm_b, nao, nao, false);
    wf_analysis_static::analysis->analyse_opdm(std::cout,
            std::string(name), std::string(desc), dm);
}


void analyse_optdm(const char *name, const char *desc,
    double *tdm_a, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;

    ab_matrix tdm(true);
    tdm.alpha() = arma::mat(tdm_a, nao, nao, false);
    wf_analysis_static::analysis->analyse_optdm(std::cout,
            std::string(name), std::string(desc), tdm);
}


void analyse_optdm(const char *name, const char *desc,
    double *tdm_a, double *tdm_b, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;

    ab_matrix tdm(false);
    tdm.alpha() = arma::mat(tdm_a, nao, nao, false);
    tdm.beta()  = arma::mat(tdm_b, nao, nao, false);
    wf_analysis_static::analysis->analyse_optdm(std::cout,
            std::string(name), std::string(desc), tdm);
}


void post_process_optdm(double *tdm_a, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;
    if (! wf_analysis_static::analysis->setup_sa_ntos(std::cout)) return;

    ab_matrix tdm(true);
    tdm.alpha() = arma::mat(tdm_a, nao, nao, false);
    wf_analysis_static::analysis->post_process_optdm(std::cout, tdm);
}


void post_process_optdm(double *tdm_a, double *tdm_b, size_t nao) {

    if (wf_analysis_static::analysis.get() == 0) return;
    if (! wf_analysis_static::analysis->setup_sa_ntos(std::cout)) return;

    ab_matrix tdm(false);
    tdm.alpha() = arma::mat(tdm_a, nao, nao, false);
    tdm.beta()  = arma::mat(tdm_b, nao, nao, false);
    wf_analysis_static::analysis->post_process_optdm(std::cout, tdm);
}


void reset_analysis() {

    if (wf_analysis_static::analysis.get() == 0) return;

    wf_analysis_static::analysis->reset();
}


void shutdown_analysis() {

    if (wf_analysis_static::analysis.get() == 0) return;

    delete wf_analysis_static::analysis.release();

}
