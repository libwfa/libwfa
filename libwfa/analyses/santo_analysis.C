#include "santo_analysis.h"
#include "nto_analysis.h"


namespace libwfa {
    
using namespace arma;

santo_analysis::santo_analysis(const arma::Mat<double> &s, const ab_matrix &c,
         const ab_matrix &edm, const ab_matrix &hdm,
         const ev_printer_i &evpr, export_data_i &opr,
         std::ostream &out) :
         m_uinv(hdm.is_alpha_eq_beta()), m_vinv_t(edm.is_alpha_eq_beta()) {

    nto_analysis_basic nab(s, c, edm, hdm);
    nab.perform(evpr, opr, out);
    
    ab_matrix ab_s(true); ab_s.alpha() = s;
    
    m_uinv   = nab.get_eigvect(false).t() * ab_s;
    m_vinv_t = ab_s * nab.get_eigvect(true);
}

void santo_analysis::stave_nto_header(std::ostream& out)
{
    out << "  Printing decomposition for individual states:\n"
        << "   H-0, H-1, ... occupied (hole) NTOs\n"
        << "   L+0, L+1, ... virtual (particle) NTOs\n"
        << " Note: only contributions due to single excitations are considered."
        << std::endl;
}

void santo_analysis::print(const ab_matrix& x, std::ostream& out)
{
    if (x.is_alpha_eq_beta()) {
        print(x.alpha(), out, 2.);
    }
    else {
        out << " Alpha:" << std::endl;
        print(x.alpha(), out);
    
        out << " Beta:" << std::endl;
        print(x.beta(), out);
    }
}

void santo_analysis::print(const arma::Mat< double >& xm, std::ostream& out,
    double dfac)
{
    char str_buf[100];


    double prt_thr = 1e-2;
    
    vec  xmv = vectorise(xm % xm);
    uvec indices = sort_index(xmv, "descend");
    
    vec col0 = xm.col(0);
    size_t nrows = col0.size();
    
    for (int i = 0; i < 100; i++) {
        size_t imat = indices(i);
        
        int imo = nrows-1 - imat % nrows;
        int jmo = nrows-1 - imat / nrows;
        double coeff = xm(imat);
        double wt = coeff * coeff * dfac;
        
        if (wt < prt_thr) break;
        
        sprintf(str_buf, "    H-%2i -> L+%2i: % .4f (%4.1f%%)\n",imo,jmo,coeff,wt*100);
        out << str_buf;
    }
        sprintf(str_buf, "                   omega = %4.1f%%\n",
            accu(xm%xm) * 100 * dfac);
    out << str_buf;

}

    
} // namespace libwfa




