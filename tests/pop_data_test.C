#include <cstdlib>
#include <sstream>
#include <libwfa/libwfa_exception.h>
#include <libwfa/analyses/pop_data.h>
#include "pop_data_test.h"

namespace libwfa {


void pop_data_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_exc();
}


void pop_data_test::test_1() {

    // Basic functionality test

    static const char *testname = "pop_data_test::test_1()";

    size_t na = 4;
    std::vector<std::string> labels(na, "");
    labels[0] = "C";
    labels[1] = "H";
    labels[2] = "H";
    labels[3] = "H";

    pop_data p;
    p.add("Set 1").randu(na);
    p.add("Set 2").randu(na);
    p.add("Set 3").randu(na);

    std::stringstream ss;
    p.print(ss, labels);

    // Check for proper line length
    std::string line;
    while (std::getline(ss, line)) {
        if (line.length() > 80) {
            fail_test(testname, __FILE__, __LINE__, "Line to long.");
        }
    }
}


void pop_data_test::test_2() {

    // Test adjustment of columns

    static const char *testname = "pop_data_test::test_2()";

    size_t na = 4;
    std::vector<std::string> labels(na, "");
    labels[0] = "C";
    labels[1] = "H";
    labels[2] = "H";
    labels[3] = "H";

    pop_data p;
    p.add("Set 1").randu(na);
    p.add("Set 2").randu(na);
    p.add("Set 3").randu(na);
    p.add("Set 4").randu(na);
    p.add("Set 5").randu(na);

    std::stringstream ss;
    p.print(ss, labels);

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


void pop_data_test::test_exc() {

    // Test exception

    static const char *testname = "pop_data_test::test_exc()";

    size_t na = 4;
    std::vector<std::string> labels(na, "");
    labels[0] = "C";
    labels[1] = "H";
    labels[2] = "H";
    labels[3] = "H";

    pop_data p1, p2, p3;
    p1.add("Set 1").randu(na + 1);
    p2.add("Set 2").randu(na);
    p2.add("Set 3").randu(na + 1);
    p3.add("Set 4").randu(na);
    p3.add("Set 5").randu(na - 1);

    std::ostringstream oss;

    bool ok = false;
    try {
        p1.print(oss, labels);
    } catch(libwfa_exception &e){
        ok = true;
    }
    if (! ok) {
        fail_test(testname, __FILE__, __LINE__, "Length (p1)");
    }
    ok = false;
    try {
        p2.print(oss, labels);
    } catch(libwfa_exception &e){
        ok = true;
    }
    if (! ok) {
        fail_test(testname, __FILE__, __LINE__, "Length (p2)");
    }
    ok = false;
    try {
        p3.print(oss, labels);
    } catch(libwfa_exception &e){
        ok = true;
    }
    if (! ok) {
        fail_test(testname, __FILE__, __LINE__, "Length (p3)");
    }
}


} // namespace libwfa
