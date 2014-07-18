#include "santo_analysis.h"


namespace libwfa {

santo_analysis::santo_analysis(const arma::Mat<double> &s, const ab_matrix &c,
         const ab_matrix &edm, const ab_matrix &hdm,
         const ev_printer_i &evpr, export_data_i &opr,
         std::ostream &out) :
         m_uinv_t(c.is_alpha_eq_beta()), m_vinv(c.is_alpha_eq_beta()) {

    nto_analysis_basic nab(s, c, edm, hdm);
    nab.perform(evpr, opr, out);
    
    ab_matrix abs(true); abs.alpha() = s;
    
    m_uinv_t =  abs * nab.get_eigvect(false);
    m_vinv   = (abs * nab.get_eigvect(true)).t();
}
    
} // namespace libwfa




