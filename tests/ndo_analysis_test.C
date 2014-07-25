#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/export/ev_printer_ndo.h>
#include <libwfa/export/export_data_print.h>
#include <libwfa/analyses/ndo_analysis.h>
#include "ndo_analysis_test.h"
#include "test01_data.h"
#include "test02_data.h"


namespace libwfa {

using namespace arma;


void ndo_analysis_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
}

template<typename TestData>
void ndo_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "ndo_analysis_test::test_1()";

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

        ab_matrix dm0(data.aeqb());
        dm0.alpha() = Mat<double>(nao, nao);
        if (! data.aeqb()) dm0.beta() = Mat<double>(nao, nao);
        read_ab_matrix(data, testname, "dm0", dm0);

        ev_printer_ndo evpr;
        std::ostringstream ssdat;
        export_data_print exdat(ssdat, "Test printer");
        
        for (size_t istate = 1; istate <= data.nstates(); istate++) {
            ab_matrix dm(data.aeqb());

            dm.alpha() = Mat<double>(nao, nao);
            if (! data.aeqb()) dm.beta() = Mat<double>(nao, nao);

            std::ostringstream ssdm; ssdm << "dm" << istate;
            read_ab_matrix(data, testname, ssdm.str().c_str(), dm);
            
            std::ostringstream outdel;           
            ab_matrix att, det, u;
            ab_vector ev;
            
            ab_matrix ddm = dm0 - dm;
            ndo_analysis ndoa(s, c, ddm, evpr);
            ndoa.perform(att, det, u, ev, exdat, outdel);
            
            // Activate if print-out is required:
            //std::cout << std::endl << ssdat.str() << std::endl;

            { // test alpha
                const Mat<double> &dm_x = ddm.alpha();
                const Mat<double> &u_x = u.alpha();
                const Col<double> &ev_x = ev.alpha();
                const Mat<double> &att_x = att.alpha();
                const Mat<double> &det_x = det.alpha();
                
                if (accu(abs(u_x.t() * s * dm_x * s * u_x - diagmat(ev_x)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(u_x.t() * s * u_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(dm_x - u_x * diagmat(ev_x) * u_x.t()) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(att_x + det_x - dm_x) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");                
            }
            { // test beta
                const Mat<double> &dm_x = ddm.beta();
                const Mat<double> &u_x = u.beta();
                const Col<double> &ev_x = ev.beta();
                const Mat<double> &att_x = att.beta();
                const Mat<double> &det_x = det.beta();
                
                if (accu(abs(u_x.t() * s * dm_x * s * u_x - diagmat(ev_x)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(u_x.t() * s * u_x - eye(nmo, nmo)) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(dm_x - u_x * diagmat(ev_x) * u_x.t()) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");
                if (accu(abs(att_x + det_x - dm_x) > 1e-12) != 0)
                    fail_test(testname, __FILE__, __LINE__, "Bad transform.");                
            }

            
        }
        
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}

} // namespace libwfa
