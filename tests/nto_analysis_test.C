#include <iomanip>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/export/ev_printer_nto.h>
#include <libwfa/export/export_data_none.h>
#include "nto_analysis_test.h"
#include "test01_data.h"
#include "test02_data.h"


namespace libwfa {

using namespace arma;


void nto_analysis_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
    fail_test("nto_analysis_test::perform()", __FILE__, __LINE__, "NIY");
}

template<typename TestData>
void nto_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "nto_analysis_test::test_1()";

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

        ev_printer_nto evpr;
        export_data_none exdat;

        // main loop
        for (size_t istate = 1; istate <= data.nstates(); istate++) {
            ab_matrix tdm(data.aeqb());

            tdm.alpha() = Mat<double>(nao, nao);
            if (! data.aeqb()) tdm.beta() = Mat<double>(nao, nao);

            std::ostringstream ssdm; ssdm << "tdm" << istate;
            read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);

            ab_matrix edm, hdm;
            std::ostringstream outdel;

         // for test purposes: use the nto_analysis class to create edm and hdm
            nto_analysis(s, c, tdm, evpr).perform(edm, hdm, exdat, outdel);

         // the eigenvalues and -vectors can only be accessed from nto_analysis_basic
            nto_analysis_basic nab(s, c, edm, hdm);

            ab_matrix lamh(data.aeqb());
            lamh.alpha() = diagmat(nab.get_eigval(false).alpha());
            if (! data.aeqb()) lamh.beta() = diagmat(nab.get_eigval(false).beta());

            ab_matrix lame(data.aeqb());
            lame.alpha() = diagmat(nab.get_eigval(true).alpha());
            if (! data.aeqb()) lame.beta() = diagmat(nab.get_eigval(true).beta());

            if (accu(abs(lamh.alpha() - lame.alpha()) > 1e-14) != 0)
                fail_test(testname, __FILE__, __LINE__, "hole != electron");
            if (accu(abs(lamh.beta() - lame.beta()) > 1e-14) != 0)
                fail_test(testname, __FILE__, __LINE__, "hole != electron");

            ab_matrix &u = nab.get_eigvect(false);
            ab_matrix &v = nab.get_eigvect(true);

            ab_matrix tdm_check = u * lamh * v.t();

            std::cout << "\nstate: " << istate << " tdm" << std::endl;
            tdm.alpha().print();
            std::cout << "tdm_check:" << std::endl;
            tdm_check.alpha().print();
            std::cout << "u:" << std::endl;
            u.alpha().print();
            std::cout << "lamh" << std::endl;
            lamh.alpha().print();
            std::cout << "v:" << std::endl;
            v.alpha().print();

            if (accu(abs(tdm_check.alpha() - tdm.alpha()) > 1e-14) != 0)
                fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            if (accu(abs(tdm_check.beta() - tdm.beta()) > 1e-14) != 0)
                fail_test(testname, __FILE__, __LINE__, "Bad transform.");

            /*ab_matrix u_ref(data.aeqb());
            u_ref.alpha() = Mat<double>(nao, nmo);
            if (! data.aeqb()) u_ref.beta() = Mat<double>(nao, nmo);
            std::ostringstream ssu; ssu << "u" << istate;
            read_ab_matrix(data, testname, ssu.str().c_str(), u_ref);

            ab_matrix u_diff = nab.get_eigvect(false) - u_ref;

            for (size_t iao = 0; iao < nao; iao++)
                for (size_t imo = 0; imo < nmo; imo++) {
                    if abs(u_diff.alpha()(iao, imo)) > 1e-14 {
                        std::cout << ""
                    }
                }*/

         /*   ab_matrix lam_ref(data.aeqb());
            lam_ref.alpha() = Mat<double>(1, nmo);
            if (! data.aeqb()) lam_ref.beta() = Mat<double>(1, nmo);
            std::ostringstream sslam; sslam << "lam" << istate;
            read_ab_matrix(data, testname, sslam.str().c_str(), lam_ref);

            Mat<double> &lam_ref_a = lam_ref.alpha();
            Mat<double> &lam_e_a = nab.get_eigval(true).alpha();

            for (size_t i = 0; i < nao * nmo; i++) {
                if (fabs(lam_ref_a[i] - lam_e_a[i]) > 1e-14) {
                    std::cout << "\nSV (alpha, elec) at position " << i << ": " <<
                            lam_e_a[i] << " , ref: " << lam_ref_a[i] << std::endl;
                    fail_test(testname, __FILE__, __LINE__, "wrong SV");
                }
            }

            */

            // fptmp
            /*std::cout << "state: " << istate << " (alpha, hole)" << std::endl;
            nab.get_eigvect(false).alpha().print();
            std::cout << "state: " << istate << " (alpha, elec)" << std::endl;
            nab.get_eigvect(true).alpha().print();*/

        }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


} // namespace libwfa
