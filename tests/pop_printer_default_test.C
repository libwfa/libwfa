#include <cstdlib>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/export/pop_printer_default.h>
#include "pop_printer_default_test.h"

namespace libwfa {

namespace {

void init_random_vector(size_t n, std::vector<double> &v) {

    v.resize(n, 0.0);
    for (std::vector<double>::iterator i = v.begin(); i != v.end(); i++) {
#ifdef HAVE_DRAND48
        *i = ::drand48();
#else
        *i = double(::rand()) / double(RAND_MAX);
#endif
    }
}

} // unnamed namespace


void pop_printer_default_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_exc();
}


void pop_printer_default_test::test_1() {

    // Basic functionality test

    static const char *testname = "pop_printer_default_test::test_1()";

    size_t na = 4;
    std::vector<std::string> labels(na, "");
    labels[0] = "C";
    labels[1] = "H";
    labels[2] = "H";
    labels[3] = "H";

    pop_data p;
    init_random_vector(na, p.add("Set 1"));
    init_random_vector(na, p.add("Set 2"));
    init_random_vector(na, p.add("Set 3"));

    std::stringstream ss;
    pop_printer_default pr(labels);
    pr.perform(p, ss);

    // Check for proper line length
    std::string line;
    while (std::getline(ss, line)) {
        if (line.length() > 80) {
            fail_test(testname, __FILE__, __LINE__, "Line to long.");
        }
    }

    // Uncomment to inspect the output
//    std::string sep(80, '-');
//    std::cout << std::endl << sep << std::endl;
//    std::cout << oss.str();
//    std::cout << std::endl << sep << std::endl;
}


void pop_printer_default_test::test_2() {

    // Test adjustment of columns

    static const char *testname = "pop_printer_default_test::test_2()";

    size_t na = 4;
    std::vector<std::string> labels(na, "");
    labels[0] = "C";
    labels[1] = "H";
    labels[2] = "H";
    labels[3] = "H";

    pop_data p;
    init_random_vector(na, p.add("Set 1"));
    init_random_vector(na, p.add("Set 2"));
    init_random_vector(na, p.add("Set 3"));
    init_random_vector(na, p.add("Set 4"));
    init_random_vector(na, p.add("Set 5"));

    std::stringstream ss;
    pop_printer_default pr(labels);
    pr.perform(p, ss);

    // Check for proper line length
    std::string line;
    while (std::getline(ss, line)) {
        if (line.length() > 80) {
            fail_test(testname, __FILE__, __LINE__, "Line to long.");
        }
    }

    // Uncomment to inspect the output
    //    std::string sep(80, '-');
    //    std::cout << std::endl << sep << std::endl;
    //    std::cout << oss.str();
    //    std::cout << std::endl << sep << std::endl;
}


void pop_printer_default_test::test_exc() {

    // Test exception

    static const char *testname = "pop_printer_default_test::test_exc()";

    size_t na = 4;
    std::vector<std::string> labels(na, "");
    labels[0] = "C";
    labels[1] = "H";
    labels[2] = "H";
    labels[3] = "H";

    pop_data p1, p2, p3;
    init_random_vector(na + 1, p1.add("Set 1"));
    init_random_vector(na, p2.add("Set 2"));
    init_random_vector(na + 1, p2.add("Set 3"));
    init_random_vector(na, p3.add("Set 4"));
    init_random_vector(na - 1, p3.add("Set 5"));

    std::ostringstream oss;
    pop_printer_default pr(labels);

    bool ok = false;
    try {
        pr.perform(p1, oss);
    } catch(libwfa_exception &e){
        ok = true;
    }
    if (! ok) {
        fail_test(testname, __FILE__, __LINE__, "Length (p1)");
    }
    ok = false;
    try {
        pr.perform(p2, oss);
    } catch(libwfa_exception &e){
        ok = true;
    }
    if (! ok) {
        fail_test(testname, __FILE__, __LINE__, "Length (p2)");
    }
    ok = false;
    try {
        pr.perform(p3, oss);
    } catch(libwfa_exception &e){
        ok = true;
    }
    if (! ok) {
        fail_test(testname, __FILE__, __LINE__, "Length (p3)");
    }
}


} // namespace libwfa
