#include <libwfa/core/transformations_dm.h>
#include "transformations_dm_test.h"

namespace libwfa {


void transformations_dm_test::perform() throw(libtest::test_exception) {

    test_form_eh_1a();
    test_form_eh_1b();
    test_form_om_1a();
    test_form_om_1b();
    test_diagonalize_dm_1a();
    test_diagonalize_dm_1b();
    test_form_ad_1a();
    test_form_ad_1b();
    test_form_ad_2();
}


using namespace arma;

namespace {

// Test data
class system_data {
private:
    static const size_t k_nao = 8;
    static const size_t k_nmo = 7;
    static const double k_s[64];
    static const double k_ca[56];
    static const double k_cb[56];

public:
    static size_t dim_ao() { return k_nao; }
    static size_t dim_mo() { return k_nmo; }

    static ab_matrix mo_coeff(bool aeqb) {
        ab_matrix c(aeqb);
        c.alpha() = Mat<double>(k_ca, k_nao, k_nmo);
        if (! aeqb) c.beta() = Mat<double>(k_cb, k_nao, k_nmo);
        return c;
    }

    static Mat<double> overlap() {
        Mat<double> s(k_s, k_nao, k_nao);
        return s;
    }

};


const double system_data::k_s[64] = {
    1.0000000000000, 0.4798622026574, 0.4247240815312, 0.5125277324114,
    0.5591394654475, 0.5424749886151, 0.5638401331380, 0.0729317713995,
    0.4798622026574, 1.0000000000000, 0.6218805948738, 0.5208388050087,
    0.3683298227843, 0.1633256506175, 0.3504020639230, 0.6737392391078,
    0.4247240815312, 0.6218805948738, 1.0000000000000, 0.7444724270608,
    0.4742539422587, 0.5067293227185, 0.5579860177822, 0.6205487444531,
    0.5125277324114, 0.5208388050087, 0.7444724270608, 1.0000000000000,
    0.4881381776650, 0.5081636528485, 0.5177969194483, 0.4249032372609,
    0.5591394654475, 0.3683298227843, 0.4742539422587, 0.4881381776650,
    1.0000000000000, 0.5588486928027, 0.4658763930202, 0.2816681417171,
    0.5424749886151, 0.1633256506175, 0.5067293227185, 0.5081636528485,
    0.5588486928027, 1.0000000000000, 0.6756330339704, 0.2136637554504,
    0.5638401331380, 0.3504020639230, 0.5579860177822, 0.5177969194483,
    0.4658763930202, 0.6756330339704, 1.0000000000000, 0.3106911403593,
    0.0729317713995, 0.6737392391078, 0.6205487444531, 0.4249032372609,
    0.2816681417171, 0.2136637554504, 0.3106911403593, 1.0000000000000
};


const double system_data::k_ca[56] = {
    -1.4386374805937,  1.7273835631264, -0.2481673864552,  0.0479179178774,
     0.1356497801719,  0.8318525719125,  0.1341667273838, -1.3595739082117,
     0.0850381455299, -0.0562263732423,  1.7516616217993, -1.1503457280556,
     0.0585263117149, -0.1472533590961, -0.1813524100028, -0.6238546456024,
    -0.4760120852191, -0.1899719329269,  0.0595910518036,  0.2555389280814,
     0.6069018646424, -1.2527765306048,  1.1466886027800, -0.3323401870186,
     0.1825908555139,  0.3537458281904, -0.4163563635022, -1.0065610904879,
    -0.1297034217118,  0.1260619715565,  0.7116904877333,  0.4311809416185,
     0.3517351768970,  0.1519289305386,  0.1593021360851,  0.3399518455318,
    -1.1239729292126, -0.1325893320848,  0.4367665279563, -0.3464214785548,
    -0.7174843291853, -0.5581910006327,  0.2570407372517,  0.0572963601103,
    -0.2632793560247,  0.5000885529216,  0.3621625819052,  0.3205582354551,
    -0.2912347029489,  0.3856273338740,  0.1900512993005,  0.0458631483621,
    -0.1963887340220, -0.3586375137324, -0.2187854906981,  0.5044324883450
};


const double system_data::k_cb[56] = {
    -1.4315480636661,  1.7109043326431, -0.5485700147404,  0.2469453763977,
     0.1234259678329,  0.8447851396450,  0.1636199488363, -1.2305877032621,
     0.0941847259226,  0.3068026509167,  1.4268212119308, -1.1016570609886,
    -0.2331361624209,  0.6258976043250, -0.7078375717425, -0.5703537020889,
    -0.4983850514309, -0.0483953845576,  0.8997159057882, -0.3233522180844,
     0.5683660872675, -1.0872543845810,  0.9028516011049, -0.7204638008000,
     0.0969889158937,  0.3136542552823, -0.4250793749127, -1.0609160516551,
     0.1165027678816,  0.1323294914388,  0.6208504711875,  0.4815469922544,
     0.5529520893900,  0.3588986380182,  0.0035099544455,  0.1030540375835,
    -1.0178677481495, -0.2331960141615,  0.4650438513266, -0.3261913045147,
    -0.5350252079374, -0.5413603936305,  0.2299895114151,  0.0776460766286,
    -0.5029167531278,  0.5106056382204,  0.5299488880724,  0.1531558750939,
    -0.4226045092075,  0.2229669218637,  0.1524084386925, -0.0076240622465,
    -0.3268928895579, -0.3139279781949, -0.1809968809003,  0.4597507108957
};


} // unnamed namespace


void transformations_dm_test::test_form_eh_1a() {

    static const char *testname = "transformations_dm_test::test_form_eh_1a()";

    //
    // Test transform of restricted TDM into electron and hole densities
    //

    size_t nb = system_data::dim_ao();

    Mat<double> s = system_data::overlap();

    ab_matrix tdm(nb, nb), de, dh;
    ab_matrix de_ref(nb, nb), dh_ref(nb, nb);

    tdm.alpha().randu();
    { // Compute reference data

    const Mat<double> &tdm_a = tdm.alpha();
    Mat<double> &de_a = de_ref.alpha(), &dh_a = dh_ref.alpha();
    for (size_t i = 0; i < nb; i++)
    for (size_t j = 0; j < nb; j++) {

        double tmp1 = 0.0, tmp2 = 0.0;
        for (size_t k = 0; k < nb; k++)
        for (size_t l = 0; l < nb; l++) {
            tmp1 += tdm_a(i, k) * s(k, l) * tdm_a(j, l);
            tmp2 += tdm_a(k, i) * s(k, l) * tdm_a(l, j);
        }
        de_a(i, j) = tmp1;
        dh_a(i, j) = tmp2;
    }

    } // End computing reference data

    // Perform operation
    try {

    form_eh(s, tdm, de, dh);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (de.nrows_a() != nb || de.ncols_a() != nb) {
        fail_test(testname, __FILE__, __LINE__, "Dim(de)");
    }
    if (dh.nrows_a() != nb || dh.ncols_a() != nb) {
        fail_test(testname, __FILE__, __LINE__, "Dim(dh)");
    }
    if (accu(abs(de.alpha() - de.alpha().t()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "de non-symmetric.");
    }
    if (accu(abs(dh.alpha() - dh.alpha().t()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "dh non-symmetric.");
    }
    if (accu(abs(de.alpha() - de_ref.alpha()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "de does not match reference.");
    }
    if (accu(abs(dh.alpha() - dh_ref.alpha()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "dh does not match reference.");
    }
}


void transformations_dm_test::test_form_eh_1b() {

    static const char *testname = "transformations_dm_test::test_form_eh_1b()";

    //
    // Test transform of unrestricted TDM into electron and hole densities
    //

    size_t nb = system_data::dim_ao();

    Mat<double> s = system_data::overlap();

    ab_matrix tdm(nb, nb, nb), de, dh;
    ab_matrix de_ref(nb, nb, nb), dh_ref(nb, nb, nb);

    tdm.alpha().randu();
    tdm.beta().randu();
    { // Compute reference data

    const Mat<double> &tdm_a = tdm.alpha(), &tdm_b = tdm.beta();
    Mat<double> &de_a = de_ref.alpha(), &de_b = de_ref.beta();
    Mat<double> &dh_a = dh_ref.alpha(), &dh_b = dh_ref.beta();
    for (size_t i = 0; i < nb; i++)
    for (size_t j = 0; j < nb; j++) {

        double tmp1a = 0.0, tmp1b = 0.0, tmp2a = 0.0, tmp2b = 0.0;
        for (size_t k = 0; k < nb; k++)
        for (size_t l = 0; l < nb; l++) {
            tmp1a += tdm_a(i, k) * s(k, l) * tdm_a(j, l);
            tmp1b += tdm_b(i, k) * s(k, l) * tdm_b(j, l);
            tmp2a += tdm_a(k, i) * s(k, l) * tdm_a(l, j);
            tmp2b += tdm_b(k, i) * s(k, l) * tdm_b(l, j);
        }
        de_a(i, j) = tmp1a;
        de_b(i, j) = tmp1b;
        dh_a(i, j) = tmp2a;
        dh_b(i, j) = tmp2b;
    }

    } // End computing reference data

    // Perform operation
    try {

    form_eh(s, tdm, de, dh);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (de.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "de: alpha == beta");
    }
    if (dh.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "dh: alpha == beta");
    }
    if (accu(abs(de.alpha() - de.alpha().t()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "de(alpha) non-symmetric.");
    }
    if (accu(abs(de.beta() - de.beta().t()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "de(beta) non-symmetric.");
    }
    if (accu(abs(dh.alpha() - dh.alpha().t()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "dh(alpha) non-symmetric.");
    }
    if (accu(abs(dh.beta() - dh.beta().t()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "dh(beta) non-symmetric.");
    }
    if (accu(abs(de.alpha() - de_ref.alpha()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "de(alpha) does not match reference.");
    }
    if (accu(abs(de.beta() - de_ref.beta()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "de(beta) does not match reference.");
    }
    if (accu(abs(dh.alpha() - dh_ref.alpha()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "dh(alpha) does not match reference.");
    }
    if (accu(abs(dh.beta() - dh_ref.beta()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "dh(beta) does not match reference.");
    }
}


void transformations_dm_test::test_form_om_1a() {

    static const char *testname = "transformations_dm_test::test_form_om_1a()";

    //
    // Test transform of restricted TDM into CT number matrix
    //

    size_t nb = system_data::dim_ao();

    Mat<double> s = system_data::overlap();

    ab_matrix tdm(nb, nb), om, om_ref(nb, nb);

    tdm.alpha().randu();
    { // Compute reference data

    const Mat<double> &tdm_a = tdm.alpha();
    Mat<double> &om_a = om_ref.alpha();
    for (size_t i = 0; i < nb; i++)
    for (size_t j = 0; j < nb; j++) {

        double tmp1 = 0.0, tmp2 = 0.0;
        for (size_t k = 0; k < nb; k++) {
            tmp1 += tdm_a(i, k) * s(k, j);
            tmp2 += s(i, k) * tdm_a(k, j);
        }
        om_a(i, j) = tmp1 * tmp2;
    }

    } // End computing reference data

    // Perform operation
    try {

    form_om(s, tdm, om);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (om.nrows_a() != nb || om.ncols_a() != nb) {
        fail_test(testname, __FILE__, __LINE__, "Dim(om)");
    }
    if (accu(abs(om.alpha() - om_ref.alpha()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "om does not match reference.");
    }
}


void transformations_dm_test::test_form_om_1b() {

    static const char *testname = "transformations_dm_test::test_form_om_1b()";

    //
    // Test transform of unrestricted TDM into electron and hole densities
    //

    size_t nb = system_data::dim_ao();

    Mat<double> s = system_data::overlap();

    ab_matrix tdm(nb, nb, nb), om, om_ref(nb, nb, nb);

    tdm.alpha().randu();
    tdm.beta().randu();
    { // Compute reference data

    const Mat<double> &tdm_a = tdm.alpha(), &tdm_b = tdm.beta();
    Mat<double> &om_a = om_ref.alpha(), &om_b = om_ref.beta();
    for (size_t i = 0; i < nb; i++)
    for (size_t j = 0; j < nb; j++) {

        double tmp1a = 0.0, tmp1b = 0.0, tmp2a = 0.0, tmp2b = 0.0;
        for (size_t k = 0; k < nb; k++) {
            tmp1a += tdm_a(i, k) * s(k, j);
            tmp1b += tdm_b(i, k) * s(k, j);
            tmp2a += s(i, k) * tdm_a(k, j);
            tmp2b += s(i, k) * tdm_b(k, j);
        }
        om_a(i, j) = tmp1a * tmp2a;
        om_b(i, j) = tmp1b * tmp2b;
    }

    } // End computing reference data

    // Perform operation
    try {

    form_om(s, tdm, om);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (om.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "om: alpha == beta");
    }
    if (accu(abs(om.alpha() - om_ref.alpha()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "om(alpha) does not match reference.");
    }
    if (accu(abs(om.beta() - om_ref.beta()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "om(beta) does not match reference.");
    }
}


void transformations_dm_test::test_diagonalize_dm_1a() {

    static const char *testname =
            "transformations_dm_test::test_diagonalize_dm_1a()";

    // Test diagonalization of a restricted and symmetric density matrix

    size_t nb = system_data::dim_ao(), nmo = system_data::dim_mo();

    Mat<double> s = system_data::overlap();

    ab_matrix dm(nb, nb), c = system_data::mo_coeff(true);

    ab_vector ev;
    ab_matrix u;

    dm.alpha().randu();
    dm.alpha() = 0.5 * (dm.alpha() + dm.alpha().t());

    // Perform operation
    try {

        diagonalize_dm(c, dm, ev, u);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (! ev.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "ev: alpha != beta");
    }
    if (! u.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "u: alpha != beta");
    }

    const Mat<double> &u_a = u.alpha(), &dm_a = dm.alpha();
    const Col<double> &ev_a = ev.alpha();
    if (accu(abs(u_a.t() * dm_a * u_a - diagmat(ev_a)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_a.t() * s * u_a - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
}


void transformations_dm_test::test_diagonalize_dm_1b() {

    static const char *testname =
            "transformations_dm_test::test_diagonalize_dm_1b()";

    // Test diagonalization of a restricted and symmetric density matrix

    size_t nb = system_data::dim_ao(), nmo = system_data::dim_mo();

    Mat<double> s = system_data::overlap();

    ab_matrix dm(nb, nb, nb), c = system_data::mo_coeff(false);

    ab_vector ev;
    ab_matrix u;

    dm.alpha().randu();
    dm.beta().randu();
    dm.alpha() = 0.5 * (dm.alpha() + dm.alpha().t());
    dm.beta()  = 0.5 * (dm.beta()  + dm.beta().t());

    // Perform operation
    try {

        diagonalize_dm(c, dm, ev, u);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (ev.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "ev: alpha == beta");
    }
    if (u.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "u: alpha == beta");
    }

    const Mat<double> &u_a = u.alpha(), &u_b = u.beta();
    const Mat<double> &dm_a = dm.alpha(), &dm_b = dm.beta();
    const Col<double> &ev_a = ev.alpha(), &ev_b = ev.beta();
    if (accu(abs(u_a.t() * dm_a * u_a - diagmat(ev_a)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_b.t() * dm_b * u_b - diagmat(ev_b)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_a.t() * s * u_a - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_b.t() * s * u_b - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
}


void transformations_dm_test::test_form_ad_1a() {

    static const char *testname = "transformations_dm_test::test_form_ad_1a()";

    // Test formation of a/d densities from eigenvectors and eigenvalues of
    // restricted and symmetric density matrix

    // Inputs
    ab_vector ev(true);
    ab_matrix u(true), c = system_data::mo_coeff(true);

    // Outputs
    ab_matrix da, dd;

    try { // Preparation of inputs

        size_t nb = system_data::dim_ao();

        Mat<double> s = system_data::overlap();
        ab_matrix dm(nb, nb);

        dm.alpha().randu();
        dm.alpha() = 0.5 * (dm.alpha() + dm.alpha().t());

        diagonalize_dm(c, dm, ev, u);
        u.alpha() = u.alpha().t() * s; // Form the transform NO2MO

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Perform operation
    try {

        form_ad(ev, u, da, dd);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (! da.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "da: alpha != beta");
    }
    if (! dd.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "dd: alpha != beta");
    }


    // Form the transform AO2NO
    Mat<double> u_a = c.alpha() * c.alpha().t() * u.alpha().t();

    Mat<double> evp;

    const Mat<double> &da_a = da.alpha(), &dd_a = dd.alpha();
    evp = u_a.t() * da_a * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density.");
    }
    evp = u_a.t() * dd_a * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density.");
    }
    evp = u_a.t() * da_a * u_a - evp;
    if (accu(abs(evp.diag() - ev.alpha()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(da) - ev(dd) != ev.");
    }
}


void transformations_dm_test::test_form_ad_1b() {

    static const char *testname = "transformations_dm_test::test_form_ad_1b()";

    // Test formation of a/d densities from unrestricted and symmetric
    // density matrix

    // Inputs
    ab_vector ev(false);
    ab_matrix u(false), c = system_data::mo_coeff(false);

    // Outputs
    ab_matrix da, dd;

    try { // Preparation of inputs

        size_t nb = system_data::dim_ao();

        Mat<double> s = system_data::overlap();
        ab_matrix dm(nb, nb, nb);

        dm.alpha().randu();
        dm.alpha() = 0.5 * (dm.alpha() + dm.alpha().t());
        dm.beta().randu();
        dm.beta() = 0.5 * (dm.beta() + dm.beta().t());

        diagonalize_dm(c, dm, ev, u);

        // Form the transform NO2MO
        u.alpha() = u.alpha().t() * s;
        u.beta() = u.beta().t() * s;

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Perform operation
    try {

        form_ad(ev, u, da, dd);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (da.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "da: alpha == beta");
    }
    if (dd.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "dd: alpha == beta");
    }

    // Form the transform AO2NO
    Mat<double> u_a = c.alpha() * c.alpha().t() * u.alpha().t();
    Mat<double> u_b = c.beta()  * c.beta().t()  * u.beta().t();

    Mat<double> evp;

    const Mat<double> &da_a = da.alpha(), &da_b = da.beta();
    evp = u_a.t() * da_a * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density (alpha).");
    }
    evp = u_b.t() * da_b * u_b;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density (beta).");
    }

    const Mat<double> &dd_a = dd.alpha(), &dd_b = dd.beta();
    evp = u_a.t() * dd_a * u_a;
    if (accu(evp % (abs(evp) > 1e-11) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density (alpha).");
    }
    evp = u_b.t() * dd_b * u_b;
    if (accu(evp % (abs(evp) > 1e-11) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density (beta).");
    }

    evp = u_a.t() * da_a * u_a - u_a.t() * dd_a * u_a;
    if (accu(abs(evp.diag() - ev.alpha()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(da) - ev(dd) != ev (alpha).");
    }
    evp = u_b.t() * da_b * u_b - u_b.t() * dd_b * u_b;
    if (accu(abs(evp.diag() - ev.beta()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(da) - ev(dd) != ev (beta).");
    }
}


void transformations_dm_test::test_form_ad_2() {

    static const char *testname = "transformations_dm_test::test_form_ad_2()";

    // Test formation of a/d densities from unrestricted and symmetric
    // density matrix

    size_t nb = system_data::dim_ao();

    Mat<double> s = system_data::overlap();

    ab_matrix dm(nb, nb, nb), c = system_data::mo_coeff(true);
    ab_matrix da, dd;

    dm.alpha().randu();
    dm.beta().randu();
    dm.alpha() = 0.5 * (dm.alpha() + dm.alpha().t());
    dm.beta() = 0.5 * (dm.beta() + dm.beta().t());

    // Perform operation
    try {

        form_ad(s, c, dm, da, dd);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Check result
    if (da.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "da: alpha == beta");
    }
    if (dd.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "dd: alpha == beta");
    }

    const Mat<double> &da_a = da.alpha(), &da_b = da.beta();
    const Mat<double> &dd_a = dd.alpha(), &dd_b = dd.beta();

    if (accu(abs(da_a - dd_a - dm.alpha()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "da - dd != dm (alpha).");
    }
    if (accu(abs(da_b - dd_b - dm.beta()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "da - dd != dm (beta).");
    }
}


} // namespace libwfa
