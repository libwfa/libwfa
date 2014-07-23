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
    //fail_test("nto_analysis_test::perform()", __FILE__, __LINE__, "NIY");
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

            ab_vector &lamh = nab.get_eigval(false);
            ab_vector &lame = nab.get_eigval(true);

            if (accu(abs(lamh.alpha() - lame.alpha()) > 1e-14) != 0)
                fail_test(testname, __FILE__, __LINE__, "hole != electron");
            if (accu(abs(lamh.beta() - lame.beta()) > 1e-14) != 0)
                fail_test(testname, __FILE__, __LINE__, "hole != electron");

            { // test alpha
                // Note: tdm apparently has to be transposed for this test to work
                const Mat<double> &tdm_x = tdm.alpha().t();
                const Mat<double> &u_x = nab.get_eigvect(false).alpha();
                const Mat<double> &v_x = nab.get_eigvect(true).alpha();
                const Col<double> &ev_x = lamh.alpha();
                Mat<double> ev_chk_x = u_x.t() * s * tdm_x * s * v_x;
                
                if (accu(abs(u_x.t() * s * u_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "U not unitary.");
                if (accu(abs(v_x.t() * s * v_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "V not unitary.");
                if (accu(abs(ev_chk_x % ev_chk_x - diagmat(ev_x)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            }
            { // test beta
                // Note: tdm apparently has to be transposed for this test to work
                const Mat<double> &tdm_x = tdm.beta().t();
                const Mat<double> &u_x = nab.get_eigvect(false).beta();
                const Mat<double> &v_x = nab.get_eigvect(true).beta();
                const Col<double> &ev_x = lamh.beta();
                Mat<double> ev_chk_x = u_x.t() * s * tdm_x * s * v_x;
                
                if (accu(abs(u_x.t() * s * u_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "U not unitary.");
                if (accu(abs(v_x.t() * s * v_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "V not unitary.");
                if (accu(abs(ev_chk_x % ev_chk_x - diagmat(ev_x)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            }
            
        }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


} // namespace libwfa
