#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/analyses/ndo_analysis.h>
#include "compare_ref.h"
#include "ndo_analysis_test.h"
#include "test_orbital_printer.h"
#include "test00_data.h"
#include "test01_data.h"
#include "test02_data.h"


namespace libwfa {

using namespace arma;


void ndo_analysis_test::perform() throw(libtest::test_exception) {

    test_form_ad_1a();
    test_form_ad_1b();
    test_form_ad_2<test01_data>();
    test_form_ad_2<test02_data>();
    test_1<test01_data>();
    test_1<test02_data>();
}



void ndo_analysis_test::test_form_ad_1a() throw(libtest::test_exception) {

    static const char *testname = "ndo_analysis_test::test_form_ad_1a()";

    // Test formation of a/d densities from eigenvectors and eigenvalues of
    // restricted and symmetric density matrix

    try { // Preparation of inputs

    size_t nao = test00_data::k_nao;
    size_t nmo = test00_data::k_nmo;
    test00_data data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(nao, nmo), dm(true), at, de;
    read_matrix(data, testname, "c_a", c.alpha());

    dm.alpha() = randu(nmo, nmo);
    dm.alpha() = c.alpha() * (dm.alpha() + dm.alpha().t()) * c.alpha().t();

    ndo_analysis(s, c, dm).form_ad(at, de);

    // Check result
    if (! at.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "da: alpha != beta");
    }
    if (! de.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "dd: alpha != beta");
    }

    orbital_data orb_a(s, c.alpha(), dm.alpha());
    const mat &u_a = orb_a.get_coeff();
    const mat &at_a = at.alpha(), &de_a = de.alpha();

    mat evp;
    evp = u_a.t() * s * at_a * s * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density.");
    }
    evp = u_a.t() * s * de_a * s * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) > 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density.");
    }
    evp = u_a.t() * s * at_a * s * u_a + evp;
    if (accu(abs(evp.diag() - orb_a.get_occ()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(da) + ev(dd) != ev.");
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


void ndo_analysis_test::test_form_ad_1b() throw(libtest::test_exception) {

    static const char *testname = "ndo_analysis_test::test_form_ad_1b()";

    // Test formation of a/d densities from unrestricted and symmetric
    // density matrix

    try { // Preparation of inputs

    size_t nao = test00_data::k_nao;
    size_t nmo = test00_data::k_nmo;
    test00_data data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(nao, nmo, nao, nmo), dm(false), at, de;
    read_matrix(data, testname, "c_a", c.alpha());
    read_matrix(data, testname, "c_b", c.beta());

    dm.alpha() = randu(nmo, nmo);
    dm.alpha() = c.alpha() * (dm.alpha() + dm.alpha().t()) * c.alpha().t();
    dm.beta() = randu(nmo, nmo);
    dm.beta() = c.beta() * (dm.beta() + dm.beta().t()) * c.beta().t();

    ndo_analysis(s, c, dm).form_ad(at, de);

    // Check result
    if (at.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "at: alpha == beta");
    }
    if (de.is_alpha_eq_beta()) {
        fail_test(testname, __FILE__, __LINE__, "de: alpha == beta");
    }

    orbital_data orb_a(s, c.alpha(), dm.alpha());
    orbital_data orb_b(s, c.beta(), dm.beta());
    const mat &u_a = orb_a.get_coeff(), &u_b = orb_b.get_coeff();
    const mat &at_a = at.alpha(), &at_b = at.beta();

    mat evp;
    evp = u_a.t() * s * at_a * s * u_a;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density (alpha).");
    }
    evp = u_b.t() * s * at_b * s * u_b;
    if (accu((evp % (abs(evp) > 1e-11)) < 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Attachment density (beta).");
    }

    const mat &de_a = de.alpha(), &de_b = de.beta();

    evp = u_a.t() * s * de_a * s * u_a;
    if (accu(evp % (abs(evp) > 1e-11) > 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density (alpha).");
    }
    evp = u_b.t() * s * de_b * s * u_b;
    if (accu(evp % (abs(evp) > 1e-11) > 0.0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Detachment density (beta).");
    }

    evp = u_a.t() * s * at_a * s * u_a + u_a.t() * s * de_a * s * u_a;
    if (accu(abs(evp.diag() - orb_a.get_occ()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(at) + ev(de) != ev (alpha).");
    }
    evp = u_b.t() * s * at_b * s * u_b + u_b.t() * s * de_b * s * u_b;
    if (accu(abs(evp.diag() - orb_b.get_occ()) > 1e-11) != 0) {
        fail_test(testname, __FILE__, __LINE__, "ev(at) + ev(de) != ev (beta).");
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


template<typename TestData>
void ndo_analysis_test::test_form_ad_2() throw(libtest::test_exception) {

    static const char *testname = "transformations_dm_test::test_form_ad_2()";

    try {

    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(data.aeqb());
    c.alpha() = mat(nao, nmo);
    if (! data.aeqb()) c.beta() = mat(nao, nmo);
    read_ab_matrix(data, testname, "c", c);

    for (size_t i = 1; i <= data.nstates(); i++) {

        ab_matrix ddm(data.aeqb());
        ab_matrix at_ref(data.aeqb()), de_ref(data.aeqb());
        ddm.alpha() = mat(nao, nao);
        at_ref.alpha() = mat(nao, nao);
        de_ref.alpha() = mat(nao, nao);
        if (! data.aeqb()) {
            ddm.beta() = mat(nao, nao);
            at_ref.beta() = mat(nao, nao);
            de_ref.beta() = mat(nao, nao);
        }
        read_ab_matrix(data, testname, "ddm1", ddm);
        read_ab_matrix(data, testname, "atdm1", at_ref);
        read_ab_matrix(data, testname, "dedm1", de_ref);

        ab_matrix at, de;
        ndo_analysis(s, c, ddm).form_ad(at, de);

        compare_ref::compare(testname, at, at_ref, 1e-14);
        compare_ref::compare(testname, de, de_ref, 1e-14);
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


/** \brief Tests the generation of attachment / detachment densities
 **/
template<typename TestData>
void ndo_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "ndo_analysis_test::test_1()";

    try {

    double thresh = 1e-1;
    test_orbital_printer pr(thresh, orbital_type::flag_t(orbital_type::NDO));

    std::ofstream of("ndo_analysis_test", std::ofstream::app);
    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;


    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(data.aeqb());
    c.alpha() = mat(nao, nmo);
    if (! data.aeqb()) c.beta() = mat(nao, nmo);
    read_ab_matrix(data, testname, "c", c);

    ab_matrix dm0(data.aeqb());
    dm0.alpha() = mat(nao, nao);
    if (! data.aeqb()) dm0.beta() = mat(nao, nao);
    read_ab_matrix(data, testname, "dm0", dm0);

    for (size_t istate = 1; istate <= data.nstates(); istate++) {

        ab_matrix dm(data.aeqb()), ddm(data.aeqb());
        dm.alpha() = mat(nao, nao);
        if (! data.aeqb()) dm.beta() = mat(nao, nao);

        std::ostringstream ssdm; ssdm << "dm" << istate;
        read_ab_matrix(data, testname, ssdm.str().c_str(), dm);

        ddm.alpha() = dm0.alpha() - dm.alpha();
        if (! data.aeqb()) ddm.beta() = dm0.beta() - dm.beta();

        ab_matrix at, de;
        ndo_analysis na(s, c, ddm);
        na.form_ad(at, de);

        const mat &ddm_a = ddm.alpha();
        const mat &at_a = at.alpha();
        const mat &de_a = de.alpha();
        orbital_data orb_a(s, c.alpha(), ddm_a);
        const mat &u_a = orb_a.get_coeff();
        const vec &e_a = orb_a.get_occ();

        if (accu(abs(u_a.t() * s * ddm_a * s * u_a - diagmat(e_a)) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform.");
        if (accu(abs(u_a.t() * s * u_a - eye(nmo, nmo)) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform.");
        if (accu(abs(ddm_a - u_a * diagmat(e_a) * u_a.t()) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform.");
        if (accu(abs(at_a + de_a - ddm_a) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform.");

        if (ddm.is_alpha_eq_beta()) {

            const mat &ddm_b = ddm.beta();
            const mat &at_b = at.beta();
            const mat &de_b = de.beta();
            orbital_data orb_b(s, c.alpha(), ddm_a);
            const mat &u_b = orb_b.get_coeff();
            const vec &e_b = orb_b.get_occ();

            if (accu(abs(u_b.t() * s * ddm_b * s * u_b - diagmat(e_b)) > 1e-12) != 0)
                fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            if (accu(abs(u_b.t() * s * u_b - eye(nmo, nmo)) > 1e-12) != 0)
                fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            if (accu(abs(ddm_b - u_b * diagmat(e_b) * u_b.t()) > 1e-12) != 0)
                fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            if (accu(abs(at_b + de_b - ddm_b) > 1e-12) != 0)
                fail_test(testname, __FILE__, __LINE__, "Bad transform.");
        }

        na.analyse(of);
        of << std::endl;
        na.export_orbitals(pr, thresh);
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}

} // namespace libwfa
