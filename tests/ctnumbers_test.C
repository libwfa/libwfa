#include <iomanip>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include "ctnumbers_test.h"
#include "test01_data.h"
#include "test02_data.h"
//fptmp
#include <libwfa/export/ctnum_export.h>

namespace libwfa {

using namespace arma;


void ctnumbers_test::perform() throw(libtest::test_exception) {
    test_1<test01_data>();
    test_1<test02_data>();
}

template<typename TestData>
void ctnumbers_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "ctnumbers_test::test_1()";

    try {
        
        size_t nao = TestData::k_nao;
        size_t nmo = TestData::k_nmo;
        size_t na  = TestData::k_natoms;
        TestData data;
        Col<size_t> b2p = data.bf2nuclei();
        
        ctnum_analysis cta(b2p);

        Mat<double> s(nao, nao);
        read_matrix(data, testname, "s", s);

        for (size_t istate = 1; istate <= data.nstates(); istate++) {
            ab_matrix tdm(data.aeqb());
            tdm.alpha() = Mat<double>(nao, nao);
            if (! data.aeqb()) tdm.beta() = Mat<double>(nao, nao);

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
            om_at_ref.alpha() = Mat<double>(na, na);
            if (! data.aeqb()) om_at_ref.beta() = Mat<double>(na, na);
            
            std::ostringstream ssom; ssom << "om" << istate;
            read_ab_matrix(data, testname, ssom.str().c_str(), om_at_ref);
                        
            for (size_t i = 0; i < na * na; i++) {
                if (fabs(om_at.alpha()[i] - om_at_ref.alpha()[i]) > 1e-14) {
                    std::ostringstream oss;
                    oss << "\nCT number (alpha) of atom " << i / na << " and atom " << i % na <<
                        "(diff: " << std::setprecision(6) << std::scientific <<
                        om_at.alpha()[i] - om_at_ref.alpha()[i] << ")";
                fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
                }
            }
            
            if (! data.aeqb()) {
                if (om_at.beta().n_rows != na || om_at.beta().n_cols != na) {
                    fail_test(testname, __FILE__, __LINE__, "Size of CT number data");
                }
                            
                for (size_t i = 0; i < na * na; i++) {
                    if (fabs(om_at.beta()[i] - om_at_ref.beta()[i]) > 1e-14) {
                        std::ostringstream oss;
                        oss << "\nCT number (beta) of atom " << i / na << " and atom " << i % na <<
                            "(diff: " << std::setprecision(6) << std::scientific <<
                            om_at.beta()[i] - om_at_ref.beta()[i] << ")";
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
