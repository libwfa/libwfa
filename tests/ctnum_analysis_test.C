#include <iomanip>
#include <libwfa/analyses/ctnum_analysis.h>
#include "ctnum_analysis_test.h"
#include "test00_data.h"

namespace libwfa {

using namespace arma;

void ctnum_analysis_test::perform() throw(libtest::test_exception) {

    test_form_om_1();
    test_1();
}


void ctnum_analysis_test::test_form_om_1() throw(libtest::test_exception) {

    static const char *testname = "ctnum_analysis_test::test_form_om_1()";

    try {

    //
    // Test transform of restricted TDM into CT number matrix
    //

    size_t nao = test00_data::k_nao;
    test00_data data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    mat tdm(nao, nao), om, om_ref(nao, nao);

    tdm.randu();
    { // Compute reference data

    for (size_t i = 0; i < nao; i++)
    for (size_t j = 0; j < nao; j++) {

        double ds_ij = 0.0, sd_ij = 0.0;
        double sds_ij = 0.0;
        for (size_t k = 0; k < nao; k++) {
            ds_ij += tdm(i, k) * s(k, j);
            sd_ij += s(i, k) * tdm(k, j);

            double ds_kj = 0.0;
            for (size_t l = 0; l < nao; l++) {
                ds_kj += tdm(k, l) * s(l, j);
            }
            sds_ij += s(i, k) * ds_kj;
        }
        om_ref(i, j) = 0.5 * (ds_ij * sd_ij + tdm(i, j) * sds_ij);
    }

    } // End computing reference data

    // Perform operation

    ctnum_analysis::form_om(s, s, tdm, "mulliken", om);

    // Check result
    if (om.n_rows != nao || om.n_cols != nao) {
        fail_test(testname, __FILE__, __LINE__, "Dim(om)");
    }
    if (accu(abs(om - om_ref) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "om does not match reference.");
    }

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


void ctnum_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "ctnum_analysis_test::test_1()";

    try {

    test00_data data;
    size_t nao = test00_data::k_nao;
    size_t na = 2, nao1 = 6;

    std::vector<uword> b2c(nao, 0);
    for (size_t i = nao1; i < nao; i++) b2c[i] = 1;

    // Use the upper and lower triagonal of a random matrix to
    // form omega matrices
    mat base = randu< mat >(nao, nao);
    mat tdm = symmatl(base);
    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    mat om_ref(na, na), om_ao;
    om_ref.fill(0.0);
    ctnum_analysis::form_om(s, s, tdm, "mulliken", om_ao);
    for (size_t i = 0; i < nao1; i++) {
        for (size_t j = 0; j < nao1; j++) {
            om_ref(0, 0) += om_ao(i, j);
        }
    }
    for (size_t i = 0; i < nao1; i++) {
        for (size_t j = nao1; j < nao; j++) {
            om_ref(0, 1) += om_ao(i, j);
        }
    }
    for (size_t i = nao1; i < nao; i++) {
        for (size_t j = 0; j < nao1; j++) {
            om_ref(1, 0) += om_ao(i, j);
        }
    }
    for (size_t i = nao1; i < nao; i++) {
        for (size_t j = nao1; j < nao; j++) {
            om_ref(1, 1) += om_ao(i, j);
        }
    }

    {
        mat om_at;
        ctnum_analysis(s, b2c, "mulliken").perform(tdm, om_at);

        if (om_at.n_rows != na || om_at.n_cols != na) {
            fail_test(testname, __FILE__, __LINE__, "Size of CT number data (Mulliken)");
        }
        for (size_t i = 0; i < na * na; i++) {
            if (fabs(om_at[i] - om_ref[i]) > 1e-14) {
                std::ostringstream oss;
                oss << "CT number of atom " << i / na << " and atom " << i % na <<
                        "(diff: " << std::setprecision(6) << std::scientific <<
                        om_at[i] - om_ref[i] << ")";
                fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
            }
        }
    }
    {
        mat om_at;
        ctnum_analysis(s, b2c, "lowdin").perform(tdm, om_at);

        if (om_at.n_rows != na || om_at.n_cols != na) {
            fail_test(testname, __FILE__, __LINE__, "Size of CT number data (Lowdin)");
        }
    }
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
