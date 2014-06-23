#include <libwfa/analyses/ex_analyse.h>
#include <libwfa/core/multipol_con.h>
#include <libwfa/core/transformations_dm.h>
#include <libwfa/export/ex_ana_printer_ad.h>
#include <fstream>
#include "ex_ana_pri_ad_test.h"


namespace libwfa{

using namespace arma;
void ex_ana_pri_ad_test::perform() throw(libtest::test_exception){
    static const char *testname = "ex_ana_pri_ad_test::test()";

        arma::Mat<double> s(4,4,fill::randu);
        arma::Mat<double> mx(4,4,fill::randu);
        arma::Mat<double> my(4,4,fill::randu);
        arma::Mat<double> mz(4,4,fill::randu);
        arma::Mat<double> mxx(4,4,fill::randu);
        arma::Mat<double> myy(4,4,fill::randu);
        arma::Mat<double> mzz(4,4,fill::randu);

        ab_matrix det(true);
        arma::Mat<double> tmp1(4,4,fill::randu);
        det.alpha()=tmp1;

        ab_matrix att(true);
        arma::Mat<double> tmp2(4,4,fill::randu);
        att.alpha()=tmp2;

        multipol_con con(mx, mxx, my, myy, mz, mzz, s);

        ex_analyse_ad analyse;
        analyse.perform(att,det,con);

        std::ofstream file("ex_ana_pri_ad_test.txt");
        ex_ana_printer_ad pr;
        pr.perform(analyse,file);


}//end fct

}//end namespace libwfa

