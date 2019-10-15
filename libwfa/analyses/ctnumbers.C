#include <iomanip>
#include "ctnumbers.h"

namespace libwfa {


using namespace arma;


ctnumbers::ctnumbers(const ctnum_analysis_i &a, const ab_matrix &tdm) :
    m_om(tdm.is_alpha_eq_beta()) {

    a.perform(tdm.alpha(), m_om.alpha());
    m_tot[0] = m_tot[1] = accu(m_om.alpha());
    m_Phe[0] = m_Phe[1] = accu(tdm.alpha() % tdm.alpha().t());
//    std::cout << "tmp: " << trace(tdm.alpha() * tdm.alpha())  << std::endl;
//    std::cout << "tmp2: " <<  accu(tdm.alpha() % tdm.alpha().t()) << std::endl;

    if (! m_om.is_alpha_eq_beta()) {
        a.perform(tdm.beta(), m_om.beta());
        m_tot[1] = accu(m_om.beta());
        m_Phe[1] = accu(tdm.beta() % tdm.beta().t());
    }
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

    out << "  <Phe> = ";
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
