#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/export/ev_printer_nto.h>
#include <libwfa/export/export_data_none.h>
#include <libwfa/core/transformations_dm.h>
#include "sa_nto_analysis_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void sa_nto_analysis_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
    //fail_test("sa_nto_analysis_test::perform()", __FILE__, __LINE__, "NIY");
}

template<typename TestData>
void sa_nto_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "sa_nto_analysis_test::test_1()";

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

        ab_matrix hdms(data.aeqb()), edms(data.aeqb());
        hdms.alpha() = mat(nao, nao, fill::zeros);
        edms.alpha() = mat(nao, nao, fill::zeros);
        if (! data.aeqb()) {
            hdms.beta() = mat(nao, nao, fill::zeros);
            edms.beta() = mat(nao, nao, fill::zeros);
        }
        
        // Preparation loop to form summed hole and particle densities
        for (size_t istate = 1; istate <= data.nstates(); istate++) {
            ab_matrix tdmt(data.aeqb());
            tdmt.alpha() = mat(nao, nao);
            if (! data.aeqb()) tdmt.beta() = mat(nao, nao);
            std::ostringstream ssdm; ssdm << "tdm" << istate;
            read_ab_matrix(data, testname, ssdm.str().c_str(), tdmt);

            ab_matrix tdm = tdmt.t();

            ab_matrix hdm, edm;
            form_eh(s, tdm, edm, hdm);
            
            hdms += hdm;
            edms += edm;
        }

        std::ostringstream outdel2;
        ev_printer_nto evpr;
        export_data_none exdat;
        sa_nto_analysis sana(s, c, edms, hdms, evpr, exdat, outdel2);
        
        // spin averaged analysis
        ab_matrix edmss(true), hdmss(true);
        edmss.alpha() = 0.5 * (edms.alpha() + edms.beta());
        hdmss.alpha() = 0.5 * (hdms.alpha() + hdms.beta());
        sa_nto_analysis ssana(s, c, edmss, hdmss, evpr, exdat, outdel2);
       
        // main loop
        for (size_t istate = 1; istate <= data.nstates(); istate++) {
            ab_matrix tdmt(data.aeqb());
            tdmt.alpha() = mat(nao, nao);
            if (! data.aeqb()) tdmt.beta() = mat(nao, nao);
            std::ostringstream ssdm; ssdm << "tdm" << istate;
            read_ab_matrix(data, testname, ssdm.str().c_str(), tdmt);

            ab_matrix tdm = tdmt.t();

            std::ostringstream outdel;
            
            { // SA-NTO decomposition for one state separately
                //std::cout << "\nsummary for state " << istate << std::endl;
                ab_matrix hdm, edm;
                form_eh(s, tdm, edm, hdm);
                sa_nto_analysis sana_i(s, c, edm, hdm, evpr, exdat, outdel);        
                
                //std::cout << "\nindividual for state " << istate << std::endl;
                ab_matrix x;
                sana_i.decompose(tdm, x);
                sana_i.print(x, outdel);
                
                check(tdm, sana_i.get_trans(false), sana_i.get_trans(true),
                      x, s, testname);
                
                for (int irow = 0; irow < x.nrows_a(); irow++) {
                    for (int jcol = 0; jcol < x.ncols_a(); jcol++) {
                        if (jcol != irow) {
                            double x_a_ij = x.alpha()(irow, jcol);
                            double x_b_ij = x.beta()(irow, jcol);
                            if (x_a_ij * x_a_ij + x_b_ij * x_b_ij > 1e-20) {
                                std::cout << "not diagonal: " << irow << " " <<
                                    jcol << std::endl;
                                //fail_test(testname, __FILE__, __LINE__, "Not diagonal.");
                            }
                        }
                    }
                }
            }
            
            { // SA-NTO decomposition for the state-averaged matrices
                //std::cout << "\nstate-ave. for state " << istate << std::endl;
                ab_matrix x;
                sana.decompose(tdm, x);
                sana.print(x, outdel);
                
                check(tdm, sana.get_trans(false), sana.get_trans(true),
                      x, s, testname);
            }
            
            { // SA-NTO decomposition for spin- and state-averaged matrices
                //std::cout << "\nspin- and state-ave. for state " << istate << std::endl;
                ab_matrix x;
                ssana.decompose(tdm, x);
                ssana.print(x, outdel);
                
                ab_matrix &ui = ssana.get_trans(false);
                ab_matrix &vit = ssana.get_trans(true);
                check(tdm, ui, vit, x, s, testname);
                
                if ( (! ui.is_alpha_eq_beta()) || (! vit.is_alpha_eq_beta()) ) {
                    fail_test(testname, __FILE__, __LINE__, "alpha != beta");
                }
            }
        }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


void sa_nto_analysis_test::check(const ab_matrix &tdm,
    const ab_matrix &ui, const ab_matrix &vit, const ab_matrix &x,
    const mat &s, const char* testname) {

    check(tdm.alpha(), ui.alpha(), vit.alpha(), x.alpha(), s, testname);
    check(tdm.beta(),  ui.beta(),  vit.beta() , x.beta(), s, testname);
}

void sa_nto_analysis_test::check(const mat &tdm_x, const mat &ui_x, const mat &vit_x,
           const mat &x_x, const mat &s, const char* testname) {

    mat x_chk_x = ui_x * tdm_x * vit_x;

    if (accu(abs(ui_x.t() * ui_x - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "U not unitary.");
    if (accu(abs(vit_x * vit_x.t() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "V not unitary.");
    if (accu(abs(x_chk_x - x_x) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "Bad transform.");

}


} // namespace libwfa
