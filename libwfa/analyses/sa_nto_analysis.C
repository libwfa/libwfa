#include <iomanip>
#include "sa_nto_analysis.h"


namespace libwfa {

using namespace arma;

sa_nto_analysis::sa_nto_analysis(const arma::mat &s, const nto_analysis &nto) :
     m_ul(nto.is_alpha_eq_beta()), m_ur(nto.is_alpha_eq_beta()) {

    if (nto.is_alpha_eq_beta()) {
        const orbital_data &nto_e = nto.get_ntos(true, false);
        const orbital_data &nto_h = nto.get_ntos(false, false);

        m_ul.alpha() = nto_h.get_coeff().t() * s;
        m_ur.alpha() = s * nto_e.get_coeff();
    }
    else {
        const orbital_data &nto_ea = nto.get_ntos(true,  false);
        const orbital_data &nto_ha = nto.get_ntos(false, false);
        const orbital_data &nto_eb = nto.get_ntos(true,  true);
        const orbital_data &nto_hb = nto.get_ntos(false, true);

        m_ul.alpha() = nto_ha.get_coeff().t() * s;
        m_ur.alpha() = s * nto_ea.get_coeff();
        m_ul.beta() = nto_hb.get_coeff().t() * s;
        m_ur.beta() = s * nto_eb.get_coeff();
    }
}


void sa_nto_analysis::decompose(const ab_matrix& tdm, ab_matrix& xdm) const {

    if (tdm.is_alpha_eq_beta()) {
        xdm.set_alpha_eq_beta();
        xdm.alpha() = m_ul.alpha() * tdm.alpha() * m_ur.alpha();
    }
    else {
        xdm.set_alpha_neq_beta();
        xdm.alpha() = m_ul.alpha() * tdm.alpha() * m_ur.alpha();
        xdm.beta() = m_ul.beta()  * tdm.beta()  * m_ur.beta();
    }
}


void sa_nto_analysis::analyse(std::ostream& out, const ab_matrix& tdm,
    ab_matrix &xdm, double thresh) const {

    decompose(tdm, xdm);
    if (tdm.is_alpha_eq_beta()) {
        analysis(out, xdm.alpha(), 2.0, 1e-2);
    }
    else {
        out << "Alpha:" << std::endl;
        analysis(out, xdm.alpha(), 1.0, 1e-2);

        out << "Beta:" << std::endl;
        analysis(out, xdm.beta(),  1.0, 1e-2);
    }
}


void sa_nto_analysis::analysis(std::ostream& out, const mat& x, double c,
    double thresh) {

    double w = accu(x%x) * c;

    out << std::fixed;
    if (w > thresh) {
        uvec indices = sort_index(vectorise(x % x), "descend");

        std::string offset(2, ' ');
        for (uword i = 0; i < indices.size(); i++) {

            uword imat = indices(i);
            uword imo = x.n_rows - 1 - imat % x.n_rows;
            uword jmo = x.n_rows - 1 - imat / x.n_rows;

            double coeff = x(imat);
            double wt = coeff * coeff * c;

            if (wt < thresh) break;

            out << offset;
            out << "H-" << std::setw(2) << imo
                    << " -> L+" << std::setw(2) << jmo << ": ";
            out << std::setw(7) << std::setprecision(4) << coeff << " ("
                    << std::setw(5) << std::setprecision(1) << wt * 100. << "%)"
                    << std::endl;
        }
    }
    out << std::string(17, ' ') << "omega = "
            << std::setw(5) << std::setprecision(1) << w * 100.
            << "%" << std::endl;
}

    
} // namespace libwfa




