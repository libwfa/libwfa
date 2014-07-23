#ifndef EXCITON_PRINTER_AD_TEST_H
#define EXCITON_PRINTER_AD_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {

/** \brief Tests for libwfa::exciton_printer_ad class

    \ingroup libwfa_tests
 **/
class exciton_printer_ad_test: public libtest::unit_test {
public:
    /** \brief Tests printer_ad and saves output in exciton_printer_ad.txt
     **/
    virtual void perform() throw(libtest::test_exception);
private:

};


} // namespace libwfa

#endif // EXCITON_PRINTER_TEST_AD_H
