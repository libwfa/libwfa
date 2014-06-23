#include <libwfa/analyses/ex_analyse.h>
#include <libwfa/core/multipol_con.h>
#include <libwfa/core/transformations_dm.h>
#include <libwfa/export/ex_ana_printer.h>
#include <fstream>
#include "ex_ana_pri_test.h"


namespace libwfa{

using namespace arma;
void ex_ana_pri_test::perform() throw(libtest::test_exception){
    static const char *testname = "ex_ana_pri_test::test()";

        arma::Mat<double> s(4,4,fill::randu);

        ab_matrix tdm(true);
        arma::Mat<double> tmp1(4,4,fill::randu);
        tdm.alpha()=tmp1;
        arma::Mat<double> mx(4,4,fill::randu);
        arma::Mat<double> my(4,4,fill::randu);
        arma::Mat<double> mz(4,4,fill::randu);
        arma::Mat<double> mxx(4,4,fill::randu);
        arma::Mat<double> myy(4,4,fill::randu);
        arma::Mat<double> mzz(4,4,fill::randu);

        ab_matrix om(true);
        arma::Mat<double> tmp2(4,4,fill::randu);
        om.alpha()=tmp2;

        multipol_con con(mx, mxx, my, myy, mz, mzz, s);

        ex_analyse analyse;
        analyse.perform(tdm, om, con);

        std::ofstream file("ex_ana_pri_test.txt");
        ex_ana_printer pr;
        pr.perform(analyse,file);

}//end fct

}//end namespace libwfa
