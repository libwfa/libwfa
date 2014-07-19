#include "transformations_dm.h"

namespace libwfa {

using namespace arma;

void form_eh(const Mat<double> &s, const ab_matrix &tdm,
    ab_matrix &de, ab_matrix &dh) {

    ab_matrix abs(true); abs.alpha() = s;
    
    de = tdm * abs * tdm.t();
    dh = tdm.t() * abs * tdm;
}


void form_om(const Mat<double> &s,
    const ab_matrix &tdm, ab_matrix &om) {
    
    ab_matrix abs(true); abs.alpha() = s;
    
    // FP: use the new formula from JCP(2014)
    om = ((tdm * abs) % (abs * tdm)) + (tdm % (abs * tdm * abs));
    om *= 0.5;
}


void diagonalize_dm(const arma::Mat<double> &s, const ab_matrix &c,
    const ab_matrix &dm, ab_vector &ev, ab_matrix &u) {

    /*  c_a is used here only for orthogonalization.
     *   Therefore it is used for both the alpha and beta density matrices.
     */
    
    mat c_a = c.alpha();
    mat cinv_a = s * c.alpha();
   
    {
        mat evec;   
        u.alpha() = cinv_a.t() * dm.alpha() * cinv_a;
        eig_sym(ev.alpha(), evec, u.alpha());
        u.alpha() = c_a * evec;
    }
    
    if (dm.is_alpha_eq_beta()) {
        ev.set_alpha_eq_beta();
        u.set_alpha_eq_beta();
    }
    else {
        mat evec;
        u.beta() = cinv_a.t() * dm.beta() * cinv_a;
        eig_sym(ev.beta(), evec, u.beta());
        u.beta() = c_a * evec;
    }
}


void form_ad(const ab_vector &ev, const ab_matrix &u,
    ab_matrix &da, ab_matrix &dd) {

    const Col<double> &ev_a = ev.alpha();
    const Mat<double> &u_a = u.alpha();
    Mat<double> &da_a = da.alpha(), &dd_a = dd.alpha();
    Col<uword> ix = find(ev_a > 0.0, 1);

    if (ix.n_rows != 0) {
        if (ix(0) != 0) {
            Mat<double> ux;
            ux = u_a.cols(ix(0), ev_a.n_rows - 1);
            da_a = ux * diagmat(ev_a.rows(ix(0), ev_a.n_rows - 1)) * ux.t();
            ux = u_a.cols(0, ix(0) - 1);
            dd_a = ux * diagmat(ev_a.rows(0, ix(0) - 1)) * ux.t();
        }
        else {
            da_a = u_a * diagmat(ev_a) * u_a.t();
            dd_a = Mat<double>(u_a.n_cols, u_a.n_cols, fill::zeros);
        }
    }
    else {
        da_a = Mat<double>(u_a.n_cols, u_a.n_cols, fill::zeros);
        dd_a = u_a * diagmat(ev_a) * u_a.t();
    }

    if (! u.is_alpha_eq_beta()) {
        da.set_alpha_neq_beta();
        dd.set_alpha_neq_beta();

        const Col<double> &ev_b = ev.beta();
        const Mat<double> &u_b = u.beta();
        Mat<double> &da_b = da.beta(), &dd_b = dd.beta();
        Col<uword> ix = find(ev_b > 0.0, 1);

        if (ix.n_rows != 0) {
            if (ix(0) != 0) {
                Mat<double> ux;
                ux = u_b.cols(ix(0), ev_b.n_rows - 1);
                da_b = ux * diagmat(ev_b.rows(ix(0), ev_b.n_rows - 1)) * ux.t();
                ux = u_b.cols(0, ix(0) - 1);
//                dd_b = ux * diagmat(ev_b.rows(0, ix(0) - 1) * -1.) * ux.t();
                dd_b = ux * diagmat(ev_b.rows(0, ix(0) - 1)) * ux.t();
            }
            else {
                da_b = u_b * diagmat(ev_b) * u_b.t();
                dd_b = Mat<double>(u_a.n_cols, u_a.n_cols, fill::zeros);
            }
        }
        else {
            da_b = Mat<double>(u_a.n_cols, u_a.n_cols, fill::zeros);
//            dd_b = u_b * diagmat(ev_b * -1.) * u_b.t();
            dd_b = u_b * diagmat(ev_b) * u_b.t();
        }
    }
    else {
        da.set_alpha_eq_beta();
        dd.set_alpha_eq_beta();
    }
}


void form_ad(const Mat<double> &s, const ab_matrix &c,
    const ab_matrix &dm, ab_matrix &da, ab_matrix &dd) {

    ab_matrix u;
    ab_vector ev;

    diagonalize_dm(s, c, dm, ev, u);
    form_ad(ev, u, da, dd);
}

void gen_transform(const ab_matrix &u, const ab_matrix &v,
        const ab_matrix &dm, ab_matrix &x) {

    x = u * dm * v.t();
}


} // namespace libwfa




