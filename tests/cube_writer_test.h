#ifndef LIBWFA_CUBE_WRITER_TEST_H
#define LIBWFA_CUBE_WRITER_TEST_H

#include <libtest/unit_test.h>

namespace libwfa {


/** \brief Tests the libwfa::grid3d class

    \ingroup libwfa_tests
 **/
class cube_writer_test : public libtest::unit_test {
public:
    virtual void perform() throw(libtest::test_exception);

private:
    void test_1();
};


} // namespace libwfa

#endif // LIBWFA_CUBE_WRITER_TEST_H
