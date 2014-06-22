#include "transformations_dm.h"

namespace libwfa {

using namespace arma;

void form_eh(const Mat<double> &s, const ab_matrix &tdm,
    ab_matrix &de, ab_matrix &dh) {

    Mat<double> &de_a = de.alpha(), &dh_a = dh.alpha();
    const Mat<double> &td_a = tdm.alpha();
    de_a = td_a * s * td_a.t();
    dh_a = td_a.t() * s * td_a;

    if (! tdm.is_alpha_eq_beta()) {
        de.set_alpha_neq_beta();
        dh.set_alpha_neq_beta();

        Mat<double> &de_b = de.beta(), &dh_b = dh.beta();
        const Mat<double> &td_b = tdm.beta();
        de_b = td_b * s * td_b.t();
        dh_b = td_b.t() * s * td_b;
    }
    else {
        de.set_alpha_eq_beta();
        dh.set_alpha_eq_beta();
    }
}


void form_om(const Mat<double> &s,
    const ab_matrix &tdm, ab_matrix &om) {

    Mat<double> &om_a = om.alpha();
    const Mat<double> &td_a = tdm.alpha();
    om_a = (td_a * s) % (s* td_a);

    if (! tdm.is_alpha_eq_beta()) {
        om.set_alpha_neq_beta();

        Mat<double> &om_b = om.beta();
        const Mat<double> &td_b = tdm.beta();
        om_b = (td_b * s) % (s* td_b);
    }
    else {
        om.set_alpha_eq_beta();
    }
}


void diagonalize_dm(const arma::Mat<double> &s, const ab_matrix &c,
    const ab_matrix &dm, ab_vector &ev, ab_matrix &u) {

    const Mat<double> &c_a = c.alpha(), &dm_a = dm.alpha();
    Col<double> &ev_a = ev.alpha();
    Mat<double> &u_a = u.alpha();
    Mat<double> evec;

    u_a = s * c_a;
    eig_sym(ev_a, evec, u_a.t() * dm_a * u_a);
    u_a = c_a * evec;

    if (! dm.is_alpha_eq_beta()) {
        ev.set_alpha_neq_beta();
        u.set_alpha_neq_beta();

        const Mat<double> &c_b = c.beta(), &dm_b = dm.beta();
        Col<double> &ev_b = ev.beta();
        Mat<double> &u_b = u.beta();

        u_b = s * c_b;
        eig_sym(ev_b, evec, u_b.t() * dm_b * u_b);
        u_b = c_b * evec;
    }
    else {
        ev.set_alpha_eq_beta();
        u.set_alpha_eq_beta();
    }
}


void form_ad(const ab_vector &ev, const ab_matrix &u,
    ab_matrix &da, ab_matrix &dd) {

    const Col<double> &ev_a = ev.alpha();
    const Mat<double> &u_a = u.alpha();
    Mat<double> &da_a = da.alpha(), &dd_a = dd.alpha();
    Col<uword> ix = find(ev_a > 0.0, 1);

    {
        Mat<double> ux;
        ux = u_a.cols(ix(0), ev_a.n_rows - 1);
        da_a = ux * diagmat(ev_a.rows(ix(0), ev_a.n_rows - 1)) * ux.t();
        ux = u_a.cols(0, ix(0) - 1);
        dd_a = ux * diagmat(ev_a.rows(0, ix(0) - 1) * -1.) * ux.t();
    }

    if (! u.is_alpha_eq_beta()) {
        da.set_alpha_neq_beta();
        dd.set_alpha_neq_beta();

        const Col<double> &ev_b = ev.beta();
        const Mat<double> &u_b = u.beta();
        Mat<double> &da_b = da.beta(), &dd_b = dd.beta();
        Col<uword> ix = find(ev_b > 0.0, 1);

        {
            Mat<double> ux;
            ux = u_b.cols(ix(0), ev_b.n_rows - 1);
            da_b = ux * diagmat(ev_b.rows(ix(0), ev_b.n_rows - 1)) * ux.t();
            ux = u_b.cols(0, ix(0) - 1);
            dd_b = ux * diagmat(ev_b.rows(0, ix(0) - 1) * -1.) * ux.t();
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

} // namespace libwfa




