#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/export/ev_printer_nto.h>
#include <libwfa/export/export_data_none.h>
#include <libwfa/core/transformations_dm.h>
#include "santo_analysis_test.h"
#include "test01_data.h"
#include "test02_data.h"

//fptmp
#include <libwfa/analyses/nto_analysis.h>

namespace libwfa {

using namespace arma;


void santo_analysis_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    //test_1<test02_data>();
    //fail_test("santo_analysis_test::perform()", __FILE__, __LINE__, "NIY");
}

template<typename TestData>
void santo_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "santo_analysis_test::test_1()";

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

        ab_matrix hdms(data.aeqb()), edms(data.aeqb());
        hdms.alpha() = Mat<double>(nao, nao, fill::zeros);
        edms.alpha() = Mat<double>(nao, nao, fill::zeros);
        if (! data.aeqb()) {
            hdms.beta() = Mat<double>(nao, nao, fill::zeros);
            edms.beta() = Mat<double>(nao, nao, fill::zeros);
        }
        
        // Preparation loop to form summed hole and particle densities
        for (size_t istate = 1; istate <= data.nstates(); istate++) {
            ab_matrix tdm(data.aeqb());

            tdm.alpha() = Mat<double>(nao, nao);
            if (! data.aeqb()) tdm.beta() = Mat<double>(nao, nao);

            std::ostringstream ssdm; ssdm << "tdm" << istate;
            read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);
            
            ab_matrix hdm, edm;
            form_eh(s, tdm, edm, hdm);
            
            hdms += hdm;
            edms += edm;
        }

        std::ostringstream outdel2;
        ev_printer_nto evpr;
        export_data_none exdat;
        santo_analysis sana(s, c, edms, hdms, evpr, exdat, outdel2);      
        
        // spin averaged analysis
        ab_matrix edmss(true), hdmss(true);
        edmss.alpha() = 0.5 * (edms.alpha() + edms.beta());
        hdmss.alpha() = 0.5 * (hdms.alpha() + hdms.beta());
        santo_analysis ssana(s, c, edmss, hdmss, evpr, exdat, outdel2);
       
        // main loop
        for (size_t istate = 1; istate <= data.nstates(); istate++) {
//        for (size_t istate = 1; istate <= 1; istate++) {
            ab_matrix tdm(data.aeqb());

            tdm.alpha() = Mat<double>(nao, nao);
            if (! data.aeqb()) tdm.beta() = Mat<double>(nao, nao);

            std::ostringstream ssdm; ssdm << "tdm" << istate;
            read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);            
            std::ostringstream outdel;
            
            { // SA-NTO decomposition for one state separately
                //std::cout << "\nsummary for state " << istate << std::endl;
                ab_matrix hdm, edm;
                form_eh(s, tdm, edm, hdm);
                santo_analysis sana_i(s, c, edm, hdm, evpr, exdat, outdel);        
                
                //std::cout << "\nindividual for state " << istate << std::endl;
                ab_matrix x;
                sana_i.decompose(tdm, x);
                sana_i.print(x, outdel);
                
                check(tdm, sana_i.get_trans(false), sana_i.get_trans(true),
                      x, s, testname);
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
                    /*ui.alpha().print();
                    std::cout << std::endl;
                    ui.beta().print();
                    std::cout << std::endl;
                    vit.alpha().print();
                    std::cout << std::endl;
                    vit.beta().print();*/
                    fail_test(testname, __FILE__, __LINE__, "alpha != beta");
                }
            }
        }

    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}

} // namespace libwfa
