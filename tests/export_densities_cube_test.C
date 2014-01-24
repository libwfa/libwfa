#include <list>
#include <libwfa/export_densities_cube.h>
#include "export_densities_cube_test.h"

namespace libwfa {

using namespace arma;

namespace {

class export_cube_test : public export_cube_base {
private:
    struct exports {
        cube::data_type type;
        std::vector<size_t> idx;
        size_t dim;

        exports(cube::data_type t, const std::vector<size_t> &i, size_t d) :
            type(t), idx(i), dim(d) { }
    };
    std::list<exports> m_lst;

public:
    typedef std::list<exports>::const_iterator iterator;

public:
    export_cube_test() : export_cube_base(cube::grid3d()) { }

    virtual ~export_cube_test() { }

    virtual void perform(cube::data_type type, const std::vector<size_t> &idx,
            const Mat<double> &data) {

        m_lst.push_back(exports(type, idx, data.n_rows));
    }

    size_t n_exports() const { return m_lst.size(); }

    iterator begin() const { return m_lst.begin(); }

    iterator end() const { return m_lst.end(); }

    cube::data_type type(iterator i) const { return i->type; }

    const std::vector<size_t> &indexes(iterator i) const { return i->idx; }

    size_t dim(iterator i) const { return i->dim; }
};

} // unnamed namespace

void export_densities_cube_test::perform() throw(libtest::test_exception) {

    test_1();
    test_2();
    test_3();
}


void export_densities_cube_test::test_1() {

    static const char *testname = "export_densities_cube_test::test_1()";

    //
    // Test exporting a single restricted density matrix
    //

    export_cube_test core;
    export_densities_cube dm_export(core);

    ab_matrix dm(5, 5);
    dm_export.perform(dm_type::sdm, 0, dm);

    if (core.n_exports() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.type(i) != cube::data_type::sdm_a) {
        fail_test(testname, __FILE__, __LINE__, "Data type.");
    }
    if (core.indexes(i).size() != 1) {
        fail_test(testname, __FILE__, __LINE__, "# density matrix in export.");
    }
    if (core.indexes(i).at(0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Index of density matrix.");
    }
    if (core.dim(i) != 25) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix.");
    }
}


void export_densities_cube_test::test_2() {

    static const char *testname = "export_densities_cube_test::test_2()";

    //
    // Test exporting a single unrestricted density matrix
    //

    export_cube_test core;
    export_densities_cube dm_export(core);

    ab_matrix dm(5, 5, 5, 4);
    dm_export.perform(dm_type::tdm, 1, dm);

    if (core.n_exports() != 2) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.type(i) != cube::data_type::tdm_a) {
        fail_test(testname, __FILE__, __LINE__, "Data type (1).");
    }
    if (core.indexes(i).size() != 1) {
        fail_test(testname, __FILE__, __LINE__,
                "# density matrix in export (1).");
    }
    if (core.indexes(i).at(0) != 1) {
        fail_test(testname, __FILE__, __LINE__, "Index of density matrix (1).");
    }
    if (core.dim(i) != 25) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix (1).");
    }
    i++;
    if (core.type(i) != cube::data_type::tdm_b) {
        fail_test(testname, __FILE__, __LINE__, "Data type (2).");
    }
    if (core.indexes(i).size() != 1) {
        fail_test(testname, __FILE__, __LINE__,
                "# density matrix in export (2).");
    }
    if (core.indexes(i).at(0) != 1) {
        fail_test(testname, __FILE__, __LINE__, "Index of density matrix (2).");
    }
    if (core.dim(i) != 20) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix (2).");
    }
}


void export_densities_cube_test::test_3() {

    static const char *testname = "export_densities_cube_test::test_3()";

    //
    // Test exporting a single restricted density matrix
    //

    export_cube_test core;
    export_densities_cube dm_export(core);

    dm_list dms(dm_type::adm);
    dms.add(0, ab_matrix(5, 5));
    dms.add(2, ab_matrix(5, 5, 5, 4));
    dm_export.perform(dms);

    if (core.n_exports() != 2) {
        fail_test(testname, __FILE__, __LINE__, "# exports.");
    }
    export_cube_test::iterator i = core.begin();
    if (core.type(i) != cube::data_type::adm_a) {
        fail_test(testname, __FILE__, __LINE__, "Data type (1).");
    }
    if (core.indexes(i).size() != 2) {
        fail_test(testname, __FILE__, __LINE__,
                "# density matrix in export (1).");
    }
    if (core.indexes(i).at(0) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Index of density matrix (1,1).");
    }
    if (core.indexes(i).at(1) != 2) {
        fail_test(testname, __FILE__, __LINE__, "Index of density matrix (1,2).");
    }
    if (core.dim(i) != 25) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix (1).");
    }
    i++;
    if (core.type(i) != cube::data_type::adm_b) {
        fail_test(testname, __FILE__, __LINE__, "Data type (2).");
    }
    if (core.indexes(i).size() != 1) {
        fail_test(testname, __FILE__, __LINE__,
                "# density matrix in export (2).");
    }
    if (core.indexes(i).at(0) != 2) {
        fail_test(testname, __FILE__, __LINE__, "Index of density matrix (2).");
    }
    if (core.dim(i) != 20) {
        fail_test(testname, __FILE__, __LINE__, "Size of density matrix (2).");
    }
}


} // namespace libwfa
