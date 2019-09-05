#include <iomanip>
#include "dyson_analysis.h"

namespace libwfa {

using namespace arma;


void dyson_analysis::analyse(std::ostream &out, size_t ndo) const {

    // For now, do nothing here
    // TODO: Further analyses like, e. g.
    //       * print leading MO contributions
    //       * angular distributions
    return;
}


void dyson_analysis::export_orbitals(orbital_printer_i &pr, double thresh) const {

    orbital_selector s;
    build_selector(m_odata.get_occ(), thresh, s);
    pr.perform(orbital_type::dyson, m_odata, s);
}


void dyson_analysis::build_selector(const vec &e, double thresh,
    orbital_selector &sel) {

    size_t ntot = e.size();

    if (sel.n_indexes() != ntot) sel = orbital_selector(ntot);

    for (size_t i = 0; i < ntot; ++i) {
        // select all Dyson orbitals whose occupation number > thresh
        if (e(i) > thresh) {
            sel.select(false, i);
        }
    }
}

} // namespace libwfa
