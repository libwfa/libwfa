#include <iomanip>
#include "exciton_analysis_base.h"
#include <libwfa/libwfa_exception.h>

namespace libwfa {

using namespace arma;


exciton_analysis_base::exciton_analysis_base(bool aeqb, size_t maxmm) {

    if (aeqb) {
        m_mom[0] = new exciton_moments(maxmm);
        m_mom[1] = 0;
    }
    else {
        m_mom[0] = new exciton_moments(maxmm);
        m_mom[1] = new exciton_moments(maxmm);
    }
}


exciton_analysis_base::~exciton_analysis_base() {

    delete m_mom[0];
    if (m_mom[1]) delete m_mom[1];
}


void exciton_analysis_base::analyse(std::ostream &out, size_t off) const {

    print_header(out, off);

    std::string os(off, ' ');
    if (m_mom[1]) {
        exciton_moments total(m_mom[0]->n_max());
        combine(*m_mom[0], *m_mom[1], total);
        out << os << "Total:" << std::endl;
        analysis(out, total, off + 2);
        out << std::endl << os << "Alpha spin:" << std::endl;
        analysis(out, *m_mom[0], off + 2);
        out << std::endl << os << "Beta spin:" << std::endl;
        analysis(out, *m_mom[1], off + 2);
    }
    else {
        analysis(out, *m_mom[0], off + 2);
    }
}


void exciton_analysis_base::print(std::ostream &out,
        const vec &vec, size_t width) {

    out << "[";
    if (vec.n_rows > 0) out << std::setw(width) << vec(0);
    for (size_t i = 1; i < vec.n_rows; i++) {
        out << ", " << std::setw(width) << vec(i);
    }
    out << "]";
}

}// end namespace libwfa

