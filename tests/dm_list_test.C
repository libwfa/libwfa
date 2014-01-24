#include <libwfa/dm_list.h>
#include "dm_list_test.h"

namespace libwfa {

using namespace arma;


void dm_list_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_3();
    test_4();
}


void dm_list_test::test_1() {

    static const char *testname = "dm_list_test::test_1()";

    //
    // Test creating an empty density matrix list
    //

    dm_list lst(dm_type::sdm);

    if (lst.type() != dm_type::sdm) {
        fail_test(testname, __FILE__, __LINE__, "Density matrix type.");
    }
    if (lst.size() != 0) {
        fail_test(testname, __FILE__, __LINE__, "List size.");
    }
    if (lst.begin() != lst.end()) {
        fail_test(testname, __FILE__, __LINE__, "Begin != end.");
    }
}


void dm_list_test::test_2() {

    static const char *testname = "dm_list_test::test_2()";

    //
    // Test construction of non-empty density matrix list
    //

    dm_list lst(dm_type::tdm);

    size_t n = 10;
    ab_matrix m1(n, n, n), m2(n, n, n);

    lst.add(1, m1);
    lst.add(4, m2);

    if (lst.type() != dm_type::tdm) {
        fail_test(testname, __FILE__, __LINE__, "Density matrix type.");
    }
    if (lst.size() != 2) {
        fail_test(testname, __FILE__, __LINE__, "List size.");
    }
    if (! lst.index_exists(1)) {
        fail_test(testname, __FILE__, __LINE__, "Index 1 missing.");
    }
    if (! lst.index_exists(4)) {
        fail_test(testname, __FILE__, __LINE__, "Index 4 missing.");
    }
    dm_list::iterator i = lst.begin();
    if (lst.get_index(i) != 1) {
        fail_test(testname, __FILE__, __LINE__, "Index 1 missing.");
    }
    i++;
    if (lst.get_index(i) != 4) {
        fail_test(testname, __FILE__, __LINE__, "Index 4 missing.");
    }
    i++;
    if (i != lst.end()) {
        fail_test(testname, __FILE__, __LINE__, "List too long.");
    }
}


void dm_list_test::test_3() {

    static const char *testname = "dm_list_test::test_3()";

    //
    // Test deleting an element from the list
    //

    dm_list lst(dm_type::adm);

    size_t n = 10;
    ab_matrix m1(n), m2(n);

    lst.add(1, m1);
    lst.add(4, m2);

    lst.erase(1);

    if (lst.type() != dm_type::adm) {
        fail_test(testname, __FILE__, __LINE__, "Density matrix type.");
    }
    if (lst.size() != 1) {
        fail_test(testname, __FILE__, __LINE__, "List size.");
    }
    if (lst.index_exists(1)) {
        fail_test(testname, __FILE__, __LINE__, "Index 1 not deleted.");
    }
    if (! lst.index_exists(4)) {
        fail_test(testname, __FILE__, __LINE__, "Index 4 missing.");
    }
    dm_list::iterator i = lst.begin();
    if (lst.get_index(i) != 4) {
        fail_test(testname, __FILE__, __LINE__, "Index 4 missing.");
    }
    i++;
    if (i != lst.end()) {
        fail_test(testname, __FILE__, __LINE__, "List too long.");
    }
}


void dm_list_test::test_4() {

    static const char *testname = "dm_list_test::test_4()";

    //
    // Test deleting an element from the list
    //

    dm_list lst(dm_type::ddm);

    size_t n1 = 3, n2 = 2;
    ab_matrix m1(n1), m2(n2);

    lst.add(2, m1);
    lst.add(3, m2);

    lst.clear();

    if (lst.type() != dm_type::ddm) {
        fail_test(testname, __FILE__, __LINE__, "Density matrix type.");
    }
    if (lst.size() != 0) {
        fail_test(testname, __FILE__, __LINE__, "List size.");
    }
    if (lst.begin() != lst.end()) {
        fail_test(testname, __FILE__, __LINE__, "List not empty.");
    }
}


} // namespace libwfa
