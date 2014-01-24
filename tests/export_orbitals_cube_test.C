#include <list>
#include <libwfa/export_orbitals_cube.h>
#include "export_orbitals_cube_test.h"

namespace libwfa {

using namespace arma;

namespace {

class export_cube_test : public export_cube_base {
private:
    struct exports {
        std::string name;
        std::vector<size_t> idx;
        size_t dim;

        exports(const std::string &n, const std::vector<size_t> &i, size_t d) :
            name(n), idx(i), dim(d) { }
    };
    std::list<exports> m_lst;

public:
    typedef std::list<exports>::const_iterator iterator;

public:
    export_cube_test() : export_cube_base(grid3d()) { }

    virtual ~export_cube_test() { }

    virtual void perform(const std::string &name, const Mat<double> &mat) {

        std::vector<size_t> idx;
        m_lst.push_back(exports(name, idx, mat.n_rows));
    }

    virtual void perform(const std::string &prefix,
            const std::vector<size_t> &idx, const Mat<double> &vecs) {

        m_lst.push_back(exports(prefix, idx, vecs.n_rows));
    }


    size_t n_exports() const { return m_lst.size(); }

    iterator begin() const { return m_lst.begin(); }

    iterator end() const { return m_lst.end(); }

    const std::string &name(iterator i) const { return i->name; }

    const std::vector<size_t> &indexes(iterator i) const { return i->idx; }

    size_t dim(iterator i) const { return i->dim; }
};

} // unnamed namespace

void export_orbitals_cube_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_3();
}


void export_orbitals_cube_test::test_1() {

    static const char *testname = "export_orbitals_cube_test::test_1()";

    //
    // Test exporting full set of restricted orbital coefficients
    //

    export_cube_test core;
    export_orbitals_cube export_c(core, "test_1");

    size_t nb = 5;
    ab_matrix c(nb);
    ab_vector ene(nb);
    ab_selector s(nb);
    s.alpha().select_all();

    export_c.perform(orbital_type::mo, c, ene, s);

    if (core.n_exports() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.name(i) != "test_1_mo") {
        fail_test(testname, __FILE__, __LINE__, "Data type.");
    }
    if (core.indexes(i).size() != nb) {
        fail_test(testname, __FILE__, __LINE__, "# coeffs in export.");
    }
    for (size_t j = 0; j < nb; j++) {
        if (core.indexes(i).at(j) != j) {
            fail_test(testname, __FILE__, __LINE__, "Index of orbital.");
        }
    }
    if (core.dim(i) != nb) {
        fail_test(testname, __FILE__, __LINE__, "Size of coeff vectors.");
    }
}


void export_orbitals_cube_test::test_2() {

    static const char *testname = "export_orbitals_cube_test::test_2()";

    //
    // Test exporting a full set of unrestricted orbital coefficients
    //

    export_cube_test core;
    export_orbitals_cube export_c(core, "test_2");

    size_t nb = 5;
    ab_matrix c(nb, nb + 1, nb, nb - 1);
    ab_vector ene(nb + 1, nb - 1);
    ab_selector s(nb + 1, nb - 1);
    s.alpha().select_all();
    s.beta().select_all();
    export_c.perform(orbital_type::mo, c, ene, s);

    if (core.n_exports() != 2) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.name(i) != "test_2_mo_a") {
        fail_test(testname, __FILE__, __LINE__, "Data type (1).");
    }
    if (core.indexes(i).size() != nb + 1) {
        fail_test(testname, __FILE__, __LINE__,
                "# orbitals in export (1).");
    }
    for (size_t j = 0; j < nb + 1; j++) {
        if (core.indexes(i).at(j) != j) {
            fail_test(testname, __FILE__, __LINE__, "Index of orbital (1).");
        }
    }
    if (core.dim(i) != nb) {
        fail_test(testname, __FILE__, __LINE__, "Size of coeff vector (1).");
    }
    i++;
    if (core.name(i) != "test_2_mo_b") {
        fail_test(testname, __FILE__, __LINE__, "Data type (2).");
    }
    if (core.indexes(i).size() != nb - 1) {
        fail_test(testname, __FILE__, __LINE__,
                "# orbitals in export (2).");
    }
    for (size_t j = 0; j < nb - 1; j++) {
        if (core.indexes(i).at(j) != j) {
            fail_test(testname, __FILE__, __LINE__, "Index of orbital (2).");
        }
    }
    if (core.dim(i) != nb) {
        fail_test(testname, __FILE__, __LINE__, "Size of coeff vector (2).");
    }
}


void export_orbitals_cube_test::test_3() {

    static const char *testname = "export_orbitals_cube_test::test_3()";

    //
    // Test exporting a partial set of unrestricted orbitals
    //

    export_cube_test core;
    export_orbitals_cube export_c(core, "test_3");

    size_t nb = 5;
    ab_matrix c(nb, nb + 1, nb, nb - 1);
    ab_vector ene(nb + 1, nb - 1);
    ab_selector s(nb + 1, nb - 1);
    s.alpha().select(2, 4);
    s.beta().deselect_all();
    export_c.perform(orbital_type::nto, c, ene, s);

    if (core.n_exports() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.name(i) != "test_3_nto_a") {
        fail_test(testname, __FILE__, __LINE__, "Data type (1).");
    }
    if (core.indexes(i).size() != 3) {
        fail_test(testname, __FILE__, __LINE__,
                "# orbitals in export (1).");
    }
    for (size_t j = 0; j < 3; j++) {
        if (core.indexes(i).at(j) != j + 2) {
            fail_test(testname, __FILE__, __LINE__, "Index of orbital (1).");
        }
    }
    if (core.dim(i) != nb) {
        fail_test(testname, __FILE__, __LINE__, "Size of coeff vector (1).");
    }
}


} // namespace libwfa
