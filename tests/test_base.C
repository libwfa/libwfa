#include "test_base.h"

namespace libwfa {


void test_base::read_double(const test_data_base &b,
    const char *testname, const char *fname, double &d) {

    if (! b.read_double(testname, fname, d)) {
        std::ostringstream oss;
        oss << "Failed to open file: " << b.make_filename(fname) << ".";
        fail_test(testname, __FILE__, __LINE__, oss.str().c_str());
    }
}


void test_base::read_ab_matrix(const test_data_base &b,
    const char *testname, const char *name, ab_matrix &m) {

    std::string fn(name);
    fn += "_a";
    read_matrix(b, testname, fn.c_str(), m.alpha());

    if (! m.is_alpha_eq_beta()) {
        fn = std::string(name) + std::string("_b");
        read_matrix(b, testname, fn.c_str(), m.beta());
    }
}


} // namespace libwfa
