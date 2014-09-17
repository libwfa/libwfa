#include <iomanip>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include "ctnumbers_test.h"
#include "test00_data.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void ctnumbers_test::perform() throw(libtest::test_exception) {

    test_form_om_1a();
    test_form_om_1b();
    test_1<test01_data>();
    test_1<test02_data>();
}


void ctnumbers_test::test_form_om_1a() throw(libtest::test_exception) {

    static const char *testname = "transformations_dm_test::test_form_om_1a()";

    try {

    //
    // Test transform of restricted TDM into CT number matrix
    //

    size_t nao = test00_data::k_nao;
    test00_data data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix tdm(nao, nao), om, om_ref(nao, nao);

    tdm.alpha().randu();
    { // Compute reference data

    const mat &tdm_a = tdm.alpha();
    mat &om_a = om_ref.alpha();
    for (size_t i = 0; i < nao; i++)
    for (size_t j = 0; j < nao; j++) {

        double ds_ij = 0.0, sd_ij = 0.0;
        double sds_ij = 0.0;
        for (size_t k = 0; k < nao; k++) {
            ds_ij += tdm_a(i, k) * s(k, j);
            sd_ij += s(i, k) * tdm_a(k, j);

            double ds_kj = 0.0;
            for (size_t l = 0; l < nao; l++) {
                ds_kj += tdm_a(k, l) * s(l, j);
            }
            sds_ij += s(i, k) * ds_kj;
        }
        om_a(i, j) = 0.5 * (ds_ij * sd_ij + tdm_a(i, j) * sds_ij);
    }

    } // End computing reference data

    // Perform operation

    ctnumbers::form_om(s, tdm, om);

    // Check result
    if (om.nrows_a() != nao || om.ncols_a() != nao) {
        fail_test(testname, __FILE__, __LINE__, "Dim(om)");
    }
    if (accu(abs(om.alpha() - om_ref.alpha()) > 1e-14) != 0) {
        fail_test(testname, __FILE__, __LINE__, "om does not match reference.");
    }

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


void ctnumbers_test::test_form_om_1b() throw(libtest::test_exception) {

    static const char *testname = "transformations_dm_test::test_form_om_1b()";

    try {

    //
    // Test transform of unrestricted TDM into electron and hole densities
    //

    size_t nao = test00_data::k_nao;
    test00_data data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix tdm(nao, nao, nao), om, om_ref(nao, nao, nao);

    tdm.alpha().randu();
    tdm.beta().randu();
    { // Compute reference data

    const mat &tdm_a = tdm.alpha(), &tdm_b = tdm.beta();
    mat &om_a = om_ref.alpha(), &om_b = om_ref.beta();
    for (size_t i = 0; i < nao; i++)
    for (size_t j = 0; j < nao; j++) {

        double tmp1a = 0.0, tmp1b = 0.0, tmp2a = 0.0, tmp2b = 0.0;
        double sds_ija = 0.0, sds_ijb = 0.0;
        for (size_t k = 0; k < nao; k++) {
            tmp1a += tdm_a(i, k) * s(k, j);
            tmp1b += tdm_b(i, k) * s(k, j);
            tmp2a += s(i, k) * tdm_a(k, j);
            tmp2b += s(i, k) * tdm_b(k, j);


            double ds_kja = 0.0, ds_kjb = 0.0;
            for (size_t l = 0; l < nao; l++) {
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

    ctnumbers::form_om(s, tdm, om);

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
        std::cout << "\nom.beta:" << std::endl;
        om.beta().print();
        std::cout << "om_ref.beta:" << std::endl;
        om_ref.beta().print();
        fail_test(testname, __FILE__, __LINE__,
                "om(beta) does not match reference.");
    }

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


template<typename TestData>
void ctnumbers_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "ctnumbers_test::test_1()";

    try {
        
    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    size_t na  = TestData::k_natoms;
    TestData data;
    uvec b2p = data.bf2nuclei();

    ctnum_analysis cta(b2p);

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    for (size_t istate = 1; istate <= data.nstates(); istate++) {
        ab_matrix tdm(data.aeqb());
        tdm.alpha() = mat(nao, nao);
        if (! data.aeqb()) tdm.beta() = mat(nao, nao);

        std::ostringstream ssdm; ssdm << "tdm" << istate;
        read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);

        ctnumbers ctnum(cta, s, tdm);

        ab_matrix om_at(data.aeqb());
        double om_tot[2];
        ctnum.perform(om_at, om_tot);

        if (om_at.alpha().n_rows != na || om_at.alpha().n_cols != na) {
            fail_test(testname, __FILE__, __LINE__, "Size of CT number data");
        }

        ab_matrix om_at_ref(data.aeqb());
        om_at_ref.alpha() = mat(na, na);
        if (! data.aeqb()) om_at_ref.beta() = mat(na, na);

        std::ostringstream ssom; ssom << "om" << istate;
        read_ab_matrix(data, testname, ssom.str().c_str(), om_at_ref);

        for (size_t i = 0; i < na * na; i++) {
            if (fabs(om_at.alpha()[i] - om_at_ref.alpha()[i]) > 1e-14) {
                std::ostringstream oss;
                oss << "CT number (alpha) of atom " << i / na
                        << " and atom " << i % na << "(diff: "
                        << std::setprecision(6) << std::scientific
                        << om_at.alpha()[i] - om_at_ref.alpha()[i] << ")";
            fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
            }
        }

        if (! data.aeqb()) {
            if (om_at.beta().n_rows != na || om_at.beta().n_cols != na) {
                fail_test(testname, __FILE__, __LINE__,
                        "Size of CT number data");
            }
                        
            for (size_t i = 0; i < na * na; i++) {
                if (fabs(om_at.beta()[i] - om_at_ref.beta()[i]) > 1e-14) {
                    std::ostringstream oss;
                    oss << "CT number (beta) of atom " << i / na
                            << " and atom " << i % na << "(diff: "
                            << std::setprecision(6) << std::scientific
                            << om_at.beta()[i] - om_at_ref.beta()[i] << ")";
                fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
                }
            }
        }
        
    }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
    
}

} // namespace libwfa
