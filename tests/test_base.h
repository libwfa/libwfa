#ifndef LIBWFA_TEST_BASE_H
#define LIBWFA_TEST_BASE_H

#include <sstream>
#include <libtest/unit_test.h>
#include "test_data_base.h"

namespace libwfa {


/** \brief Base for unit tests

    \ingroup libwfa_tests
 **/
class test_base : public libtest::unit_test {
protected:
    template<typename T>
    void read_matrix(const test_data_base &b, const char *testname,
            const char *fname, arma::Mat<T> &m);

    void read_double(const test_data_base &b, const char *testname,
            const char *fname, double &d);

    void read_ab_matrix(const test_data_base &b, const char *testname,
            const char *name, ab_matrix &m);
};


template<typename T>
void test_base::read_matrix(const test_data_base &b,
    const char *testname, const char *fname, arma::Mat<T> &m) {

    if (! b.read_matrix(testname, fname, m)) {
        std::ostringstream oss;
        oss << "read_matrix for file " << b.make_filename(fname) << " failed.";
        fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
    }
}


} // namespace libwfa

#endif // LIBWFA_TEST_BASE_H
