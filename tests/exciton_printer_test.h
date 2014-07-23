#ifndef EXCITON_PRINTER_TEST_H
#define EXCITON_PRINTER_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {

/** \brief Tests for libwfa::exciton_printer class

    \ingroup libwfa_tests
 **/
class exciton_printer_test: public libtest::unit_test {
public:
    /** \brief Tests printer and saves output in exciton_printer.txt
     **/
    virtual void perform() throw(libtest::test_exception);
private:

};


} // namespace libwfa



#endif // EXCITON_PRINTER_TEST_H
