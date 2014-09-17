#include <iomanip>
#include "ctnumbers.h"

namespace libwfa {


using namespace arma;


ctnumbers::ctnumbers(const ctnum_analysis_i &a, const arma::mat &s,
    const ab_matrix &tdm) : m_om(tdm.is_alpha_eq_beta()) {

    ab_matrix om;
    form_om(s, tdm, om);

    initialize(a, om);
}


void ctnumbers::analyse(std::ostream &out) const {

    size_t w = 10, prec = 4;
    out << "  omega = ";
    out << std::setprecision(prec) << std::fixed;
    if (m_om.is_alpha_eq_beta()) {
        out << std::setw(w) << m_tot[0] + m_tot[1] << std::endl;
    }
    else {
        out << std::setw(w) << m_tot[0] + m_tot[1] << " (alpha: ";
        out << std::setw(w) << m_tot[0] << ", beta: ";
        out << std::setw(w) << m_tot[1] << ")" << std::endl;
    }
}


void ctnumbers::form_om(const mat &s, const ab_matrix &tdm, ab_matrix &om) {

    if (tdm.is_alpha_eq_beta()) {
        om.set_alpha_eq_beta();
        form_om(s, tdm.alpha(), om.alpha());
    }
    else {
        om.set_alpha_neq_beta();
        form_om(s, tdm.alpha(), om.alpha());
        form_om(s, tdm.beta(), om.beta());
    }
}


void ctnumbers::initialize(const ctnum_analysis_i &a, const ab_matrix &om) {

    m_tot[0] = accu(om.alpha());
    a.perform(om.alpha(), m_om.alpha());

    if (m_om.is_alpha_eq_beta()) {
        m_tot[1] = 0.0;
    }
    else {
        m_tot[1] = accu(om.beta());
        a.perform(om.beta(),  m_om.beta());
    }
}


} // namespace libwfa



