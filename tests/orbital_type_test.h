#ifndef LIBWFA_ORBITAL_TYPE_TEST_H
#define LIBWFA_ORBITAL_TYPE_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::orbital_type class

    \ingroup libwfa_tests
 **/
class orbital_type_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_TYPE_TEST_H
