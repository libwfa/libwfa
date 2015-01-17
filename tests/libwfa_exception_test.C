#include <libwfa/libwfa_exception.h>
#include "libwfa_exception_test.h"

namespace libwfa {


void libwfa_exception_test::perform() throw(libtest::test_exception) {

    bool thrown = false;

    try {
        throw libwfa_exception("libwfa_exception_test", "perform()",
                __FILE__, __LINE__, "Test");
    }
    catch (libwfa_exception &e) {

        thrown = true;
    }

    if (! thrown) {
        fail_test("libwfa_exception_test::perform()",
                __FILE__, __LINE__, "Not thrown.");
    }
}


} // namespace libwfa
