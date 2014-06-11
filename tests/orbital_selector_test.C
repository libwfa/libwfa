#include <sstream>
#include <libwfa/core/orbital_selector.h>
#include "orbital_selector_test.h"

namespace libwfa {


void orbital_selector_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_3();
}


void orbital_selector_test::test_1() {

    static const char *testname = "orbital_selector_test::test_1()";

    //
    // Test setting up a orbital_selector with no elements selected
    //

    size_t n = 10;

    orbital_selector s(n);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes.");
    }
    if (s.n_selected() != 0) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    if (! s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "None selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 0) {
        fail_test(testname, __FILE__, __LINE__, "# indexes.");
    }
}


void orbital_selector_test::test_2() {

    static const char *testname = "orbital_selector_test::test_2()";

    //
    // Test setting up a orbital_selector with some elements selected
    //

    size_t n = 10;

    orbital_selector s(n);
    s.select(true, n-3, n-1, 1, true);
    s.select(false, 0, 3);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 7) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.n_selected(true) != 3) {
        fail_test(testname, __FILE__, __LINE__, "# selected occ");
    }
    if (s.n_selected(false) != 4) {
        fail_test(testname, __FILE__, __LINE__, "# selected vir");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 7) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    for (size_t i = 0; i < 3; i++) {
        if (indexes[i] != n - i - 1)
            fail_test(testname, __FILE__, __LINE__, "i");
    }
    for (size_t i = 4; i < 7; i++) {
        if (indexes[i] != i - 4)
            fail_test(testname, __FILE__, __LINE__, "i");
    }

}


void orbital_selector_test::test_3() {

    static const char *testname = "orbital_selector_test::test_5()";

    //
    // Test deselecting a previously selected element
    //

    size_t n = 10;

    orbital_selector s(n);
    s.select(true, 2, 5);
    s.select(false, 6, 8);
    s.deselect(7);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 6) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.all_selected() || s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 6) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    for (size_t i = 0; i < 4; i++) {
        if (indexes[i] != i + 2)
            fail_test(testname, __FILE__, __LINE__, "indexes[i]");
    }
    if (indexes[4] != 6) {
        fail_test(testname, __FILE__, __LINE__, "indexes[4]");
    }
    if (indexes[5] != 8) {
        fail_test(testname, __FILE__, __LINE__, "indexes[5]");
    }
}


} // namespace libwfa
