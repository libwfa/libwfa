#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/analyses/no_analysis.h>
#include <libwfa/export/ev_printer_no.h>
#include <libwfa/export/export_data_none.h>
#include "no_analysis_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void no_analysis_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();    
}

template<typename TestData>
void no_analysis_test::test_1() throw(libtest::test_exception) {

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
        
        ev_printer_no evpr;
        export_data_none exdat;
               
        // main loop
        for (size_t istate = 1; istate <= data.nstates(); istate++) {
            ab_matrix dm(data.aeqb());

            dm.alpha() = Mat<double>(nao, nao);
            if (! data.aeqb()) dm.beta() = Mat<double>(nao, nao);

            std::ostringstream ssdm; ssdm << "dm" << istate;
            read_ab_matrix(data, testname, ssdm.str().c_str(), dm);
            
            std::ostringstream outdel;
            
            no_analysis noa(s, c, dm);
            noa.perform(evpr, exdat, outdel);

            { // test alpha
                const Mat<double> &dm_x = dm.alpha();
                const Mat<double> &u_x = noa.get_eigvect(false).alpha();
                const Col<double> &ev_x = noa.get_eigval(false).alpha();
                if (accu(abs(u_x.t() * s * dm_x * s * u_x - diagmat(ev_x)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(u_x.t() * s * u_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(dm_x - u_x * diagmat(ev_x) * u_x.t()) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            }
            { // test beta
                const Mat<double> &dm_x = dm.beta();
                const Mat<double> &u_x = noa.get_eigvect(false).beta();
                const Col<double> &ev_x = noa.get_eigval(false).beta();
                if (accu(abs(u_x.t() * s * dm_x * s * u_x - diagmat(ev_x)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(u_x.t() * s * u_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(dm_x - u_x * diagmat(ev_x) * u_x.t()) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            }
            { // test spin-traced
                Mat<double> dm_x = 0.5 * (dm.alpha() + dm.beta());
                const Mat<double> &u_x = noa.get_eigvect(true).alpha();
                const Col<double> &ev_x = noa.get_eigval(true).alpha();// + noa.get_eigval(true).alpha();
                if (accu(abs(u_x.t() * s * dm_x * s * u_x - diagmat(ev_x)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(u_x.t() * s * u_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(dm_x - u_x * diagmat(ev_x) * u_x.t()) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
            }


        }
        
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
    
}

} // namespace libwfa
