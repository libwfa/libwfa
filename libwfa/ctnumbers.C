#include <iomanip>
#include "ctnumbers.h"
#include "transformations_dm.h"

namespace libwfa {


using namespace arma;


void ctnumbers::perform(const ab_matrix &tdm,
        ab_matrix &om, std::vector<double> &om_tot) {

    om_tot.resize(2, 0.0);

    ab_matrix om_ao;
    form_om(m_s, tdm, om_ao);

    om_tot[0] = accu(om_ao.alpha());
    m_analysis.perform(om_ao.alpha(), om.alpha());

    if (tdm.is_alpha_eq_beta()) {
        om_tot[1] = om_tot[0];
    }
    else {
        om_tot[1] = accu(om_ao.beta());
        m_analysis.perform(om_ao.beta(), om.beta());
    }
}

void ctnumbers::perform(const ab_matrix &tdm, std::ostream &out,
    ctnum_print_i &pr) {

    std::vector<double> om_tot;
    ab_matrix om;

    perform(tdm, om, om_tot);

    size_t w = 10, prec = 4;
    out << "omega = ";
    out << std::setw(w) << std::setprecision(prec) << std::fixed;
    out << om_tot[0] + om_tot[1] << " (alpha: ";
    out << std::setw(w) << std::setprecision(prec) << std::fixed;
    out << om_tot[0] << ", beta: ";
    out << std::setw(w) << std::setprecision(prec) << std::fixed;
    out << om_tot[1] << ")";

    pr.perform(om);
}


} // namespace libwfa



