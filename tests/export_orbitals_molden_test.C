#include <list>
#include <libwfa/export_orbitals_molden.h>
#include <libwfa/libwfa_exception.h>
#include "export_orbitals_molden_test.h"

namespace libwfa {

using namespace arma;

namespace {

class molden_file_test : public export_molden_i {
private:
    struct exports {
        size_t norbs[5];

        exports(size_t no_a, size_t nv_a,
            size_t no_b, size_t nv_b, size_t nb) {
            norbs[0] = no_a; norbs[1] = nv_a;
            norbs[2] = no_b; norbs[3] = nv_b;
            norbs[4] = nb;
        }
    };
    std::list<exports> m_lst;

public:
    typedef std::list<exports>::const_iterator iterator;

public:
    molden_file_test() { }

    virtual ~molden_file_test() { }

    virtual void perform(const std::string &name, const ab_matrix &coeff,
            const ab_vector &ene, size_t nocc_a, size_t nocc_b) {

        const Mat<double> &ca = coeff.alpha(), &cb = coeff.beta();
        const Col<double> &ea = ene.alpha(), &eb = ene.beta();

        if (ca.n_cols != ea.n_rows || cb.n_cols != eb.n_rows) {
            throw libwfa_exception("molden_file_test", "perform()",
                    __FILE__, __LINE__, "Inconsistent # orbitals");
        }
        if (cb.n_cols != 0 && ca.n_rows != cb.n_rows) {
            throw libwfa_exception("molden_file_test", "perform()",
                    __FILE__, __LINE__, "Inconsistent # basis functions");
        }

        m_lst.push_back(exports(nocc_a, ca.n_cols - nocc_a,
                nocc_b, cb.n_cols - nocc_b, ca.n_rows));
    }

    size_t n_exports() const { return m_lst.size(); }

    iterator begin() const { return m_lst.begin(); }

    iterator end() const { return m_lst.end(); }

    size_t no_a(iterator i) const { return i->norbs[0]; }
    size_t nv_a(iterator i) const { return i->norbs[1]; }
    size_t no_b(iterator i) const { return i->norbs[2]; }
    size_t nv_b(iterator i) const { return i->norbs[3]; }
    size_t nb(iterator i) const { return i->norbs[4]; }
};

} // unnamed namespace

void export_orbitals_molden_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_3();
}


void export_orbitals_molden_test::test_1() {

    static const char *testname = "export_orbitals_molden_test::test_1()";

    //
    // Test exporting full set of restricted orbital coefficients
    //

    size_t nb = 5, no = 2;

    molden_file_test core;
    export_orbitals_molden export_c(core, "test_1", no, nb - no, no, nb - no);

    ab_matrix c(nb);
    ab_vector ene(nb);
    ab_selector s(nb);
    s.alpha().select_all();

    try {
        export_c.perform(orbital_type::mo, c, ene, s);
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    if (core.n_exports() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    molden_file_test::iterator i = core.begin();
    if (core.no_a(i) != core.no_b(i) || core.no_a(i) != no) {
        fail_test(testname, __FILE__, __LINE__, "# occupied.");
    }
    if (core.nv_a(i) != core.nv_b(i) || core.nv_a(i) != nb - no) {
        fail_test(testname, __FILE__, __LINE__, "# virtual.");
    }
    if (core.nb(i) != nb) {
        fail_test(testname, __FILE__, __LINE__, "# basis functions.");
    }
}


void export_orbitals_molden_test::test_2() {

    static const char *testname = "export_orbitals_molden_test::test_2()";

    //
    // Test exporting full set of unrestricted orbital coefficients
    //

    size_t nb = 5, no_a = 2, no_b = 1;

    molden_file_test core;
    export_orbitals_molden export_c(core, "test_2",
            no_a, nb - no_a, no_b, nb - no_b);

    ab_matrix c(nb, nb, nb, nb);
    ab_vector ene(nb, nb);
    ab_selector s(nb, nb);
    s.alpha().select_all();
    s.beta().select_all();

    try {
        export_c.perform(orbital_type::mo, c, ene, s);
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    if (core.n_exports() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    molden_file_test::iterator i = core.begin();
    if (core.no_a(i) != no_a || core.no_b(i) != no_b) {
        fail_test(testname, __FILE__, __LINE__, "# occupied.");
    }
    if (core.nv_a(i) != nb - no_a || core.nv_b(i) != nb - no_b) {
        fail_test(testname, __FILE__, __LINE__, "# virtual.");
    }
    if (core.nb(i) != nb) {
        fail_test(testname, __FILE__, __LINE__, "# basis functions.");
    }
}


void export_orbitals_molden_test::test_3() {

    static const char *testname = "export_orbitals_molden_test::test_3()";

    //
    // Test exporting partial set of unrestricted orbital coefficients
    //

    size_t nb = 5, no_a = 2, no_b = 1;

    molden_file_test core;
    export_orbitals_molden export_c(core, "test_3",
            no_a, nb - no_a, no_b, nb - no_b);

    ab_matrix c(nb, nb, nb, nb);
    ab_vector ene(nb, nb);
    ab_selector s(nb, nb);
    s.alpha().select(1, 3);
    s.beta().deselect_all();

    try {
        export_c.perform(orbital_type::mo, c, ene, s);
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

    if (core.n_exports() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    molden_file_test::iterator i = core.begin();
    if (core.no_a(i) != 1 || core.no_b(i) != 0) {
        fail_test(testname, __FILE__, __LINE__, "# occupied.");
    }
    if (core.nv_a(i) != 2 || core.nv_b(i) != 0) {
        fail_test(testname, __FILE__, __LINE__, "# virtual.");
    }
    if (core.nb(i) != nb) {
        fail_test(testname, __FILE__, __LINE__, "# basis functions.");
    }
}


} // namespace libwfa
