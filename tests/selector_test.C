#include <sstream>
#include <libwfa/selector.h>
#include "selector_test.h"

namespace libwfa {


void selector_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_3();
    test_4a();
    test_4b();
    test_5();
    test_6();
    test_7();
}


void selector_test::test_1() {

    static const char *testname = "selector_test::test_1()";

    //
    // Test setting up a selector with no elements selected
    //

    size_t n = 10;

    selector s(n);

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


void selector_test::test_2() {

    static const char *testname = "selector_test::test_2()";

    //
    // Test setting up a selector with all elements selected
    //

    size_t n = 10;

    selector s(n);
    s.select_all();

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != n) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (! s.all_selected()) {
        fail_test(testname, __FILE__, __LINE__, "All selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    for (size_t i = 0; i < n; i++) {
        if (indexes[i] != i)
            fail_test(testname, __FILE__, __LINE__, "i");
    }
}


void selector_test::test_3() {

    static const char *testname = "selector_test::test_3()";

    //
    // Test setting up a selector with a few elements selected
    //

    size_t n = 10;

    selector s(n);
    s.select(2);
    s.select(5);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 2) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.all_selected() || s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 2) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (indexes[0] != 2) {
        fail_test(testname, __FILE__, __LINE__, "indexes[0]");
    }
    if (indexes[1] != 5) {
        fail_test(testname, __FILE__, __LINE__, "indexes[1]");
    }
}


void selector_test::test_4a() {

    static const char *testname = "selector_test::test_4a()";

    //
    // Test setting up a selector with a range of elements selected
    //

    size_t n = 10;

    selector s(n);
    s.select(2, 4);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.all_selected() || s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (indexes[0] != 2) {
        fail_test(testname, __FILE__, __LINE__, "indexes[0]");
    }
    if (indexes[1] != 3) {
        fail_test(testname, __FILE__, __LINE__, "indexes[1]");
    }
    if (indexes[2] != 4) {
        fail_test(testname, __FILE__, __LINE__, "indexes[2]");
    }
}


void selector_test::test_4b() {

    static const char *testname = "selector_test::test_4b()";

    //
    // Test setting up a selector with a range of elements selected
    //

    size_t n = 10;

    selector s(n);
    s.select(2, 8, 2);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 4) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.all_selected() || s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 4) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (indexes[0] != 2) {
        fail_test(testname, __FILE__, __LINE__, "indexes[0]");
    }
    if (indexes[1] != 4) {
        fail_test(testname, __FILE__, __LINE__, "indexes[1]");
    }
    if (indexes[2] != 6) {
        fail_test(testname, __FILE__, __LINE__, "indexes[2]");
    }
    if (indexes[3] != 8) {
        fail_test(testname, __FILE__, __LINE__, "indexes[3]");
    }
}


void selector_test::test_5() {

    static const char *testname = "selector_test::test_5()";

    //
    // Test deselecting a previously selected element
    //

    size_t n = 10;

    selector s(n);
    s.select(2, 5);
    s.deselect(5);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.all_selected() || s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (indexes[0] != 2) {
        fail_test(testname, __FILE__, __LINE__, "indexes[0]");
    }
    if (indexes[1] != 3) {
        fail_test(testname, __FILE__, __LINE__, "indexes[1]");
    }
    if (indexes[2] != 4) {
        fail_test(testname, __FILE__, __LINE__, "indexes[2]");
    }
}


void selector_test::test_6() {

    static const char *testname = "selector_test::test_6()";

    //
    // Test deselecting a previous selected element
    //

    size_t n = 10;

    selector s(n);
    s.select(2, 5);
    s.deselect(5);

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.all_selected() || s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 3) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (indexes[0] != 2) {
        fail_test(testname, __FILE__, __LINE__, "indexes[0]");
    }
    if (indexes[1] != 3) {
        fail_test(testname, __FILE__, __LINE__, "indexes[1]");
    }
    if (indexes[2] != 4) {
        fail_test(testname, __FILE__, __LINE__, "indexes[2]");
    }
}


void selector_test::test_7() {

    static const char *testname = "selector_test::test_7()";

    //
    // Test deselecting all elements
    //

    size_t n = 10;

    selector s(n);
    s.select(2, 5);
    s.deselect_all();

    if (s.n_indexes() != n) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
    if (s.n_selected() != 0) {
        fail_test(testname, __FILE__, __LINE__, "# selected");
    }
    if (s.all_selected() || ! s.none_selected()) {
        fail_test(testname, __FILE__, __LINE__, "# selected.");
    }
    std::vector<size_t> indexes = s.get_selected();
    if (indexes.size() != 0) {
        fail_test(testname, __FILE__, __LINE__, "# indexes");
    }
}


} // namespace libwfa
