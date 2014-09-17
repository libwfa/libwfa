#ifndef LIBWFA_ORBITAL_DATA_TEST_H
#define LIBWFA_ORBITAL_DATA_TEST_H

#include "test_base.h"

namespace libwfa {


/** \brief Tests the libwfa::orbital_data class

    \ingroup libwfa_tests
 **/
class orbital_data_test : public test_base {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1a() throw(libtest::test_exception);
    void test_1b() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_DATA_TEST_H
