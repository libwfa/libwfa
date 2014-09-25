#include <iomanip>
#include "exciton_analysis_base.h"

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
        out << os << "alpha spin:" << std::endl;
        analysis(out, *m_mom[0], off + 2);
        out << os << "beta spin:" << std::endl;
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


void exciton_analysis_base::combine(const exciton_moments &a,
    const exciton_moments &b, exciton_moments &res) {

    size_t nmax = std::min(a.n_max(), b.n_max());
    if (nmax != res.n_max()) {
        throw 1;
    }

    double sa = accu(a.get(0, 0)), sb = accu(b.get(0, 0));
    if (sa + sb == 0.0) {
        throw 1;
    }

    for (size_t i = 0; i < nmax; i++) {
        for (size_t j = 0; j <= i; j++) {

            vec ma = a.get(i, j), mb = b.get(i, j);
            vec mc = (sa * ma + sb * mb) / (sa + sb);
            res.set(i, j, mc);
        }
    }
}

}// end namespace libwfa

