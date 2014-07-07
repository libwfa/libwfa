#include <iomanip>
#include <libwfa/core/transformations_dm.h>
#include "compare_ref.h"
#include "test01_data.h"
#include "test02_data.h"
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
    test_form_ad_3<test01_data>();
    test_form_ad_3<test02_data>();
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
        c.alpha() = Mat<double>(k_ca, k_nmo, k_nao).t();
        if (! aeqb) c.beta() = Mat<double>(k_cb, k_nmo, k_nao).t();
        return c;
    }

    static Mat<double> overlap() {
        Mat<double> s(k_s, k_nao, k_nao);
        return s;
    }

};


const double system_data::k_s[64] = {
    1.0000000000000, 0.1499364919728, 0.2756397025660, 0.3608786604600,
    0.4231694217306, 0.2169289010344, 0.1068759472109, 0.2323201323161,
    0.1499364919728, 1.0000000000000, 0.3046494444134, 0.2488729695324,
    0.3538992105750, 0.1415702020749, 0.3882910344983, 0.0813575608190,
    0.2756397025660, 0.3046494444134, 1.0000000000000, 0.3447649691952,
    0.2505632289685, 0.4164878517622, 0.3013916595373, 0.2038155294722,
    0.3608786604600, 0.2488729695324, 0.3447649691952, 1.0000000000000,
    0.2846212877193, 0.3250353743788, 0.2322411764180, 0.3108507213183,
    0.4231694217306, 0.3538992105750, 0.2505632289685, 0.2846212877193,
    1.0000000000000, 0.2960018882295, 0.0548056028783, 0.0696227924200,
    0.2169289010344, 0.1415702020749, 0.4164878517622, 0.3250353743788,
    0.2960018882295, 1.0000000000000, 0.3356770082610, 0.2525453858543,
    0.1068759472109, 0.3882910344983, 0.3013916595373, 0.2322411764180,
    0.0548056028783, 0.3356770082610, 1.0000000000000, 0.2133969775168,
    0.2323201323161, 0.0813575608190, 0.2038155294722, 0.3108507213183,
    0.0696227924200, 0.2525453858543, 0.2133969775168, 1.0000000000000
};


const double system_data::k_ca[56] = {
     0.2070241232447, -0.4812337921132, -0.1390835157406, -0.4925552636112,
     0.1829069314102, -0.8584554360221, -0.2611475692609,  0.1957609095120,
     0.1487236507165,  0.6044164668592, -0.8254101697727,  0.4375103971058,
     0.2165620001094,  0.2978310514672,  0.2400229385512,  0.1233647205885,
     0.0282560407939,  0.3279427192784, -0.4341689760692, -0.3330358136162,
     0.7041376090490,  0.2387212380805, -0.0938203132651, -0.2088638489827,
     0.1126813666373,  0.2175092824450,  0.6111827841331, -0.5934145150443,
     0.2056943267020, -0.4949719693035,  0.3125492241421,  0.8707179936394,
    -0.0876926155303,  0.3158952159234,  0.1213351444442,  0.2289245292677,
     0.1415090846135, -0.1768593872072, -0.7095217090289, -0.7119527654296,
     0.2810049866600, -0.2127116890329,  0.1916909150547,  0.5675942963097,
     0.1547130574029,  0.6590510400547,  0.1443645793922, -0.4451411987105,
    -0.5402420823171,  0.1647863064104,  0.1451250483719, -0.6062917221410,
     0.0906126960968,  0.4983472022308,  0.1540374913833,  0.5808013562576
};


const double system_data::k_cb[56] = {
    -0.4812337921132, -0.1390835157406, -0.4925552636112,  0.1829069314102,
    -0.8584554360221, -0.2611475692609,  0.0222295707007,  0.1487236507165,
     0.6044164668592, -0.8254101697727,  0.4375103971058,  0.2165620001094,
     0.2978310514672,  0.0017453764321,  0.1233647205885,  0.0282560407939,
     0.3279427192784, -0.4341689760692, -0.3330358136162,  0.7041376090490,
    -0.6493452139263, -0.0938203132651, -0.2088638489827,  0.1126813666373,
     0.2175092824450,  0.6111827841331, -0.5934145150443, -0.6950709610176,
    -0.4949719693035,  0.3125492241421,  0.8707179936394, -0.0876926155303,
     0.3158952159234,  0.1213351444442,  0.4853917932760,  0.1415090846135,
    -0.1768593872072, -0.7095217090289, -0.7119527654296,  0.2810049866600,
    -0.2127116890329,  0.4544227966897,  0.5675942963097,  0.1547130574029,
     0.6590510400547,  0.1443645793922, -0.4451411987105, -0.5402420823171,
     0.2552584308317,  0.1451250483719, -0.6062917221410,  0.0906126960968,
     0.4983472022308,  0.1540374913833,  0.5808013562576,  0.3886283343235
};


} // unnamed namespace


void transformations_dm_test::test_form_eh_1a() throw(libtest::test_exception) {

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


void transformations_dm_test::test_form_eh_1b() throw(libtest::test_exception) {

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


void transformations_dm_test::test_form_om_1a() throw(libtest::test_exception) {

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

        double ds_ij = 0.0, sd_ij = 0.0;        
        double sds_ij = 0.0;
        for (size_t k = 0; k < nb; k++) {
            ds_ij += tdm_a(i, k) * s(k, j);
            sd_ij += s(i, k) * tdm_a(k, j);

            double ds_kj = 0.0;
            for (size_t l = 0; l < nb; l++) {
                ds_kj += tdm_a(k, l) * s(l, j);
            }
            sds_ij += s(i, k) * ds_kj;
        }
        om_a(i, j) = 0.5 * (ds_ij * sd_ij + tdm_a(i, j) * sds_ij);
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


void transformations_dm_test::test_form_om_1b() throw(libtest::test_exception) {

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
        double sds_ija = 0.0, sds_ijb = 0.0;
        for (size_t k = 0; k < nb; k++) {
            tmp1a += tdm_a(i, k) * s(k, j);
            tmp1b += tdm_b(i, k) * s(k, j);
            tmp2a += s(i, k) * tdm_a(k, j);
            tmp2b += s(i, k) * tdm_b(k, j);

            
            double ds_kja = 0.0, ds_kjb = 0.0;
            for (size_t l = 0; l < nb; l++) {
                ds_kja += tdm_a(k, l) * s(l, j);
                ds_kjb += tdm_b(k, l) * s(l, j);
            }
            sds_ija += s(i, k) * ds_kja;
            sds_ijb += s(i, k) * ds_kjb;

        }
        om_a(i, j) = 0.5 * (tmp1a * tmp2a + tdm_a(i, j) * sds_ija);
        om_b(i, j) = 0.5 * (tmp1b * tmp2b + tdm_b(i, j) * sds_ijb);
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
        std::cout << "\nom.alpha:" << std::endl;
        om.alpha().print();
        std::cout << "om_ref.alpha:" << std::endl;
        om_ref.alpha().print();
        fail_test(testname, __FILE__, __LINE__,
                "om(alpha) does not match reference.");
    }
    if (accu(abs(om.beta() - om_ref.beta()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "om(beta) does not match reference.");
    }
}


void transformations_dm_test::test_diagonalize_dm_1a()
throw(libtest::test_exception) {

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

        diagonalize_dm(s, c, dm, ev, u);

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

    const Mat<double> &dm_a = dm.alpha();
    const Mat<double> &u_a = u.alpha();
    const Col<double> &ev_a = ev.alpha();
    if (accu(abs(u_a.t() * s * dm_a * s * u_a - diagmat(ev_a)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_a.t() * s * u_a - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
}


void transformations_dm_test::test_diagonalize_dm_1b()
throw(libtest::test_exception) {

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

        diagonalize_dm(s, c, dm, ev, u);

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

    const Mat<double> &dm_a = dm.alpha(), &dm_b = dm.beta();
    const Mat<double> &u_a = u.alpha(), &u_b = u.beta();
    const Col<double> &ev_a = ev.alpha(), &ev_b = ev.beta();
    if (accu(abs(u_a.t() * s * dm_a * s * u_a - diagmat(ev_a)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_b.t() * s * dm_b * s * u_b - diagmat(ev_b)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_a.t() * s * u_a - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
    if (accu(abs(u_b.t() * s * u_b - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");
    }
}


void transformations_dm_test::test_form_ad_1a() throw(libtest::test_exception) {

    static const char *testname = "transformations_dm_test::test_form_ad_1a()";

    // Test formation of a/d densities from eigenvectors and eigenvalues of
    // restricted and symmetric density matrix

    // Inputs
    Mat<double> s = system_data::overlap();
    ab_vector ev(true);
    ab_matrix u(true);

    try { // Preparation of inputs

        size_t nb = system_data::dim_ao();
        ab_matrix c = system_data::mo_coeff(true);

        ab_matrix dm(nb, nb);

        dm.alpha().randu();
        dm.alpha() = 0.5 * (dm.alpha() + dm.alpha().t());

        diagonalize_dm(s, c, dm, ev, u);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Outputs
    ab_matrix da, dd;

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

    const Mat<double> &u_a = u.alpha();
    const Mat<double> &da_a = da.alpha(), &dd_a = dd.alpha();

    Mat<double> evp;
    evp = u_a.t() * s * da_a * s * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density.");
    }
    evp = u_a.t() * s * dd_a * s * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density.");
    }
    evp = u_a.t() * s * da_a * s * u_a - evp;
    if (accu(abs(evp.diag() - ev.alpha()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(da) - ev(dd) != ev.");
    }
}


void transformations_dm_test::test_form_ad_1b() throw(libtest::test_exception) {

    static const char *testname = "transformations_dm_test::test_form_ad_1b()";

    // Test formation of a/d densities from unrestricted and symmetric
    // density matrix

    // Inputs
    Mat<double> s = system_data::overlap();
    ab_vector ev(false);
    ab_matrix u(false);

    try { // Preparation of inputs

        size_t nb = system_data::dim_ao();
        ab_matrix c = system_data::mo_coeff(false);

        ab_matrix dm(nb, nb, nb);

        dm.alpha().randu();
        dm.alpha() = 0.5 * (dm.alpha() + dm.alpha().t());
        dm.beta().randu();
        dm.beta() = 0.5 * (dm.beta() + dm.beta().t());

        diagonalize_dm(s, c, dm, ev, u);

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    // Outputs
    ab_matrix da, dd;

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

    const Mat<double> &u_a = u.alpha(), &u_b = u.beta();
    const Mat<double> &da_a = da.alpha(), &da_b = da.beta();

    Mat<double> evp;
    evp = u_a.t() * s * da_a * s * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density (alpha).");
    }
    evp = u_b.t() * s * da_b * s * u_b;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density (beta).");
    }

    const Mat<double> &dd_a = dd.alpha(), &dd_b = dd.beta();

    evp = u_a.t() * s * dd_a * s * u_a;
    if (accu(evp % (abs(evp) > 1e-11) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density (alpha).");
    }
    evp = u_b.t() * s * dd_b * s * u_b;
    if (accu(evp % (abs(evp) > 1e-11) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density (beta).");
    }

    evp = u_a.t() * s * da_a * s * u_a - u_a.t() * s * dd_a * s * u_a;
    if (accu(abs(evp.diag() - ev.alpha()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(da) - ev(dd) != ev (alpha).");
    }
    evp = u_b.t() * s * da_b * s * u_b - u_b.t() * s * dd_b * s * u_b;
    if (accu(abs(evp.diag() - ev.beta()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(da) - ev(dd) != ev (beta).");
    }
}


void transformations_dm_test::test_form_ad_2() throw(libtest::test_exception) {

    static const char *testname = "transformations_dm_test::test_form_ad_2()";

    // Test formation of a/d densities from unrestricted and symmetric
    // density matrix

    size_t nb = system_data::dim_ao(), nmo = system_data::dim_mo();
    Mat<double> s = system_data::overlap();
    ab_matrix c = system_data::mo_coeff(false);
    ab_matrix dm(false);

    dm.alpha() = randu< Mat<double> >(nmo, nmo);
    dm.beta() = randu< Mat<double> >(nmo, nmo);
    dm.alpha() = c.alpha() * 0.5 * (dm.alpha() + dm.alpha().t()) * c.alpha().t();
    dm.beta() = c.beta() * 0.5 * (dm.beta() + dm.beta().t()) *c.beta().t();

    ab_matrix da, dd;

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


template<typename TestData>
void transformations_dm_test::test_form_ad_3() throw(libtest::test_exception) {

    static const char *testname = "transformations_dm_test::test_form_ad_3()";

    try {

    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;

    Mat<double> s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(data.aeqb());
    c.alpha() = Mat<double>(nao, nmo);
    if (! data.aeqb()) c.beta() = Mat<double>(nao, nmo);
    read_ab_matrix(data, testname, "c", c);

    for (size_t i = 1; i <= data.nstates(); i++) {

        ab_matrix ddm(data.aeqb());
        ab_matrix at_ref(data.aeqb()), de_ref(data.aeqb());
        ddm.alpha() = Mat<double>(nao, nao);
        at_ref.alpha() = Mat<double>(nao, nao);
        de_ref.alpha() = Mat<double>(nao, nao);
        if (! data.aeqb()) {
            ddm.beta() = Mat<double>(nao, nao);
            at_ref.beta() = Mat<double>(nao, nao);
            de_ref.beta() = Mat<double>(nao, nao);
        }
        read_ab_matrix(data, testname, "ddm1", ddm);
        read_ab_matrix(data, testname, "atdm1", at_ref);
        read_ab_matrix(data, testname, "dedm1", de_ref);

        ab_matrix at, de;

        // Perform operation
        form_ad(s, c, ddm, at, de);

        compare_ref::compare(testname, at, at_ref, 1e-14);
        compare_ref::compare(testname, de, de_ref, 1e-14);
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


} // namespace libwfa
