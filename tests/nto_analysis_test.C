#include <iomanip>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include "nto_analysis_test.h"
#include "test_orbital_printer.h"
#include "test00_data.h"
#include "test01_data.h"
#include "test02_data.h"


namespace libwfa {

using namespace arma;


void nto_analysis_test::perform() throw(libtest::test_exception) {

    test_form_eh_1();
    test_form_eh_2();
    test_1<test01_data>();
    test_1<test02_data>();
}


void nto_analysis_test::test_form_eh_1() throw(libtest::test_exception) {

    static const char *testname = "nto_analysis_test::test_form_eh_1()";

    try {

    //
    // Test transform of restricted TDM into electron and hole densities
    //

    size_t nao = test00_data::k_nao;
    size_t nmo = test00_data::k_nmo;
    test00_data data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix tdm(nao, nao), de, dh;
    ab_matrix de_ref(nao, nao), dh_ref(nao, nao);

    tdm.alpha().randu();
    { // Compute reference data

    const mat &tdm_a = tdm.alpha();
    mat &de_a = de_ref.alpha(), &dh_a = dh_ref.alpha();
    for (size_t i = 0; i < nao; i++)
    for (size_t j = 0; j < nao; j++) {

        double tmp1 = 0.0, tmp2 = 0.0;
        for (size_t k = 0; k < nao; k++)
        for (size_t l = 0; l < nao; l++) {
            tmp1 += tdm_a(i, k) * s(k, l) * tdm_a(j, l);
            tmp2 += tdm_a(k, i) * s(k, l) * tdm_a(l, j);
        }
        dh_a(i, j) = tmp1;
        de_a(i, j) = tmp2;
    }

    } // End computing reference data

    // Perform operation
    nto_analysis::form_eh(s, tdm, de, dh);

    // Check result
    if (de.nrows_a() != nao || de.ncols_a() != nao) {
        fail_test(testname, __FILE__, __LINE__, "Dim(de)");
    }
    if (dh.nrows_a() != nao || dh.ncols_a() != nao) {
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

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


void nto_analysis_test::test_form_eh_2() throw(libtest::test_exception) {

    static const char *testname = "nto_analysis_test::test_form_eh_2()";

    try {

    //
    // Test transform of unrestricted TDM into electron and hole densities
    //

    size_t nao = test00_data::k_nao;
    size_t nmo = test00_data::k_nmo;
    test00_data data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix tdm(nao, nao, nao), de, dh;
    ab_matrix de_ref(nao, nao, nao), dh_ref(nao, nao, nao);

    tdm.alpha().randu();
    tdm.beta().randu();
    { // Compute reference data

    const mat &tdm_a = tdm.alpha(), &tdm_b = tdm.beta();
    mat &de_a = de_ref.alpha(), &de_b = de_ref.beta();
    mat &dh_a = dh_ref.alpha(), &dh_b = dh_ref.beta();
    for (size_t i = 0; i < nao; i++)
    for (size_t j = 0; j < nao; j++) {

        double tmp1a = 0.0, tmp1b = 0.0, tmp2a = 0.0, tmp2b = 0.0;
        for (size_t k = 0; k < nao; k++)
        for (size_t l = 0; l < nao; l++) {
            tmp1a += tdm_a(i, k) * s(k, l) * tdm_a(j, l);
            tmp1b += tdm_b(i, k) * s(k, l) * tdm_b(j, l);
            tmp2a += tdm_a(k, i) * s(k, l) * tdm_a(l, j);
            tmp2b += tdm_b(k, i) * s(k, l) * tdm_b(l, j);
        }
        dh_a(i, j) = tmp1a;
        dh_b(i, j) = tmp1b;
        de_a(i, j) = tmp2a;
        de_b(i, j) = tmp2b;
    }

    } // End computing reference data

    // Perform operation
    nto_analysis::form_eh(s, tdm, de, dh);

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

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


template<typename TestData>
void nto_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "nto_analysis_test::test_1()";

    try {

    double thresh = 1e-1;
    test_orbital_printer pr(thresh, orbital_type::flag_t(orbital_type::NTO));

    std::ofstream of("nto_analysis_test", std::ofstream::app);
    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(data.aeqb());
    c.alpha() = mat(nao, nmo);
    if (! data.aeqb()) c.beta() = mat(nao, nmo);
    read_ab_matrix(data, testname, "c", c);

    for (size_t istate = 1; istate <= data.nstates(); istate++) {

        ab_matrix tdm(data.aeqb());
        tdm.alpha() = mat(nao, nao);
        if (! data.aeqb()) tdm.beta() = mat(nao, nao);

        std::ostringstream ssdm; ssdm << "tdm" << istate;
        read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);

        ab_matrix edm, hdm;
        nto_analysis::form_eh(s, tdm, edm, hdm);

        orbital_data nto_ea(s, c.alpha(), edm.alpha());
        orbital_data nto_ha(s, c.alpha(), hdm.alpha());
        if (accu(abs(nto_ea.get_occ() - nto_ha.get_occ()) > 1e-14) != 0)
            fail_test(testname, __FILE__, __LINE__, "hole != electron");

        orbital_data nto_eb(s, c.beta(), edm.beta());
        orbital_data nto_hb(s, c.beta(), hdm.beta());
        if (accu(abs(nto_eb.get_occ() - nto_hb.get_occ()) > 1e-14) != 0)
            fail_test(testname, __FILE__, __LINE__, "hole != electron");

        // for test purposes: use the nto_analysis class to create edm and hdm
        nto_analysis na(s, c, edm, hdm);
        na.analyse(of);
        of << std::endl;
        na.export_orbitals(pr, thresh);
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
