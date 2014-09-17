#include "transformations_dm.h"

namespace libwfa {

using namespace arma;

void form_eh(const mat &s, const ab_matrix &tdm,
    ab_matrix &de, ab_matrix &dh) {

    form_eh(s, tdm.alpha(), de.alpha(), dh.alpha());
    if (tdm.is_alpha_eq_beta()) {
        de.set_alpha_eq_beta();
        dh.set_alpha_eq_beta();
    }
    else {
        form_eh(s, tdm.beta(), de.beta(), dh.beta());
    }
}

void form_eh(const arma::mat &s, const arma::mat &tdm,
        arma::mat &de, arma::mat &dh) {
    de = tdm.t() * s * tdm;
    dh = tdm * s * tdm.t();
}

void form_om(const mat &s,
    const ab_matrix &tdm, ab_matrix &om) {

    // FP: use the new formula from JCP(2014)
    form_om(s, tdm.alpha(), om.alpha());
    if (tdm.is_alpha_eq_beta())
        om.set_alpha_eq_beta();
    else
        form_om(s, tdm.beta(), om.beta());
}

void form_om(const arma::mat &s, const arma::mat &tdm,
        arma::mat &om) {

    om = (tdm * s) % (s * tdm);
    om += tdm % (s * tdm * s);
    om *= 0.5;
}

void diagonalize_dm(const arma::mat &s, const ab_matrix &c,
    const ab_matrix &dm, ab_vector &ev, ab_matrix &u) {

    /*  c_a is used here only for orthogonalization.
     *   Therefore it is used for both the alpha and beta density matrices.
     */
    
    const mat &c_a = c.alpha();
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

    const vec &ev_a = ev.alpha();
    const mat &u_a = u.alpha();
    mat &da_a = da.alpha(), &dd_a = dd.alpha();
    Col<uword> ix = find(ev_a > 0.0, 1);

    if (ix.n_rows != 0) {
        if (ix(0) != 0) {
            mat ux;
            ux = u_a.cols(ix(0), ev_a.n_rows - 1);
            da_a = ux * diagmat(ev_a.rows(ix(0), ev_a.n_rows - 1)) * ux.t();
            ux = u_a.cols(0, ix(0) - 1);
            dd_a = ux * diagmat(ev_a.rows(0, ix(0) - 1)) * ux.t();
        }
        else {
            da_a = u_a * diagmat(ev_a) * u_a.t();
            dd_a = mat(u_a.n_cols, u_a.n_cols, fill::zeros);
        }
    }
    else {
        da_a = mat(u_a.n_cols, u_a.n_cols, fill::zeros);
        dd_a = u_a * diagmat(ev_a) * u_a.t();
    }

    if (! u.is_alpha_eq_beta()) {
        da.set_alpha_neq_beta();
        dd.set_alpha_neq_beta();

        const vec &ev_b = ev.beta();
        const mat &u_b = u.beta();
        mat &da_b = da.beta(), &dd_b = dd.beta();
        Col<uword> ix = find(ev_b > 0.0, 1);

        if (ix.n_rows != 0) {
            if (ix(0) != 0) {
                mat ux;
                ux = u_b.cols(ix(0), ev_b.n_rows - 1);
                da_b = ux * diagmat(ev_b.rows(ix(0), ev_b.n_rows - 1)) * ux.t();
                ux = u_b.cols(0, ix(0) - 1);
//                dd_b = ux * diagmat(ev_b.rows(0, ix(0) - 1) * -1.) * ux.t();
                dd_b = ux * diagmat(ev_b.rows(0, ix(0) - 1)) * ux.t();
            }
            else {
                da_b = u_b * diagmat(ev_b) * u_b.t();
                dd_b = mat(u_a.n_cols, u_a.n_cols, fill::zeros);
            }
        }
        else {
            da_b = mat(u_a.n_cols, u_a.n_cols, fill::zeros);
//            dd_b = u_b * diagmat(ev_b * -1.) * u_b.t();
            dd_b = u_b * diagmat(ev_b) * u_b.t();
        }
    }
    else {
        da.set_alpha_eq_beta();
        dd.set_alpha_eq_beta();
    }
}


void form_ad(const mat &s, const ab_matrix &c,
    const ab_matrix &dm, ab_matrix &da, ab_matrix &dd) {

    ab_matrix u;
    ab_vector ev;

    diagonalize_dm(s, c, dm, ev, u);
    form_ad(ev, u, da, dd);
}

} // namespace libwfa




