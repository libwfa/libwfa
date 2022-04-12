#include <iomanip>
#include "ctnumbers.h"

namespace libwfa {


using namespace arma;


ctnumbers::ctnumbers(const ctnum_analysis_i &a, const ab_matrix &tdm) :
    m_om(tdm.is_alpha_eq_beta()) {

    a.perform(tdm.alpha(), m_om.alpha(), m_Phe[0], m_LOC[0]);
    m_tot[0] = m_tot[1] = accu(m_om.alpha());

    if (m_om.is_alpha_eq_beta()) {
        m_Phe[1] = m_Phe[0];
        m_LOC[1] = m_LOC[0];
    }
    else {
        a.perform(tdm.beta(), m_om.beta(), m_Phe[1], m_LOC[1]);
        m_tot[1] = accu(m_om.beta());
    }

    m_om_ab = a.compute_DSDS(tdm.alpha(), tdm.beta());
}


void ctnumbers::analyse(std::ostream &out) const {

    size_t w = 10, prec = 4;
    out << "  omega        = ";
    out << std::setprecision(prec) << std::fixed;
    if (m_om.is_alpha_eq_beta()) {
        out << std::setw(w) << m_tot[0] + m_tot[1] << std::endl;
    }
    else {
        out << std::setw(w) << m_tot[0] + m_tot[1] << " (alpha: ";
        out << std::setw(w) << m_tot[0] << ", beta: ";
        out << std::setw(w) << m_tot[1] << ")" << std::endl;
    }

    out << " 2<alpha|beta> = ";
    out << std::setw(w) << 2*m_om_ab << std::endl;

    out << "  LOC          = ";
    out << std::setw(w) << m_LOC[0] + m_LOC[1];
    if (m_om.is_alpha_eq_beta()) {
        out << std::endl;
    }
    else {
        out << " (alpha: ";
        out << std::setw(w) << m_LOC[0] << ", beta: ";
        out << std::setw(w) << m_LOC[1] << ")" << std::endl;
    }

    out << "  <Phe>        = ";
    out << std::setprecision(prec) << std::fixed;
    if (m_om.is_alpha_eq_beta()) {
        out << std::setw(w) << (m_Phe[0] + m_Phe[1]) / (m_tot[0] + m_tot[1]) << std::endl;
    }
    else {
        out << std::setw(w) << (m_Phe[0] + m_Phe[1]) / (m_tot[0] + m_tot[1]) << " (alpha: ";
        out << std::setw(w) << m_Phe[0] / m_tot[0] << ", beta: ";
        out << std::setw(w) << m_Phe[1] / m_tot[1]  << ")" << std::endl;
    }
}


} // namespace libwfa
