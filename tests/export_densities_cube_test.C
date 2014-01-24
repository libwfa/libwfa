#include <list>
#include <libwfa/export_densities_cube.h>
#include "export_densities_cube_test.h"

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

void export_densities_cube_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
}


void export_densities_cube_test::test_1() {

    static const char *testname = "export_densities_cube_test::test_1()";

    //
    // Test exporting a single restricted density matrix
    //

    export_cube_test core;
    export_densities_cube dm_export(core, "state_1");

    ab_matrix dm(5, 5);
    dm_export.perform(dm_type::state, dm);

    if (core.n_exports() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.name(i) != "state_1_sdm") {
        fail_test(testname, __FILE__, __LINE__, "Name.");
    }
    if (core.indexes(i).size() != 0) {
        fail_test(testname, __FILE__, __LINE__, "# density matrix in export.");
    }
    if (core.dim(i) != 5) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix.");
    }
}


void export_densities_cube_test::test_2() {

    static const char *testname = "export_densities_cube_test::test_2()";

    //
    // Test exporting a single unrestricted density matrix
    //

    export_cube_test core;
    export_densities_cube dm_export(core, "state_2");

    ab_matrix dm(5, 5, 5, 5);
    dm_export.perform(dm_type::transition, dm);

    if (core.n_exports() != 2) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.name(i) != "state_2_tdm_a") {
        fail_test(testname, __FILE__, __LINE__, "Name (1).");
    }
    if (core.indexes(i).size() != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "# density matrix in export (1).");
    }
    if (core.dim(i) != 5) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix (1).");
    }
    i++;
    if (core.name(i) != "state_2_tdm_b") {
        fail_test(testname, __FILE__, __LINE__, "Name (2).");
    }
    if (core.indexes(i).size() != 0) {
        fail_test(testname, __FILE__, __LINE__,
                "# density matrix in export (2).");
    }
    if (core.dim(i) != 5) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix (2).");
    }
}


} // namespace libwfa
