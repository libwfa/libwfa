#include <libwfa/export/exciton_printer_ad.h>
#include <sstream>
#include "exciton_printer_ad_test.h"
#include "exciton_test_data_hhe.h"


namespace libwfa{

using namespace arma;


void exciton_printer_ad_test::perform() throw(libtest::test_exception){

    static const char *testname = "exciton_printer_ad_test::test()";

    ab_exciton_moments mom(2, true);

    mom.alpha().set(1, 0, exciton_test_data_hhe::re());
    mom.alpha().set(2, 0, exciton_test_data_hhe::re2());
    mom.alpha().set(0, 1, exciton_test_data_hhe::rh());
    mom.alpha().set(0, 2, exciton_test_data_hhe::rh2());

    std::ostringstream oss;
    exciton_printer_ad().perform(mom, oss);

    // TODO: test output properly!
    // Reference values:
    //  |re - rh|_2 = 18.8971
    //  sh = sqrt(|rh2 - rh * rh|_1) = 0.86985
    //  se = sqrt(|re2 - re * re|_1) = 1.3446

    //std::cout << std::endl << std::endl << oss.str() << std::endl;
}


}//end namespace libwfa
