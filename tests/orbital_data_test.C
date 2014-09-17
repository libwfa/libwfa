#include <iomanip>
#include <libwfa/core/orbital_data.h>
#include "compare_ref.h"
#include "test00_data.h"
#include "orbital_data_test.h"

namespace libwfa {


void orbital_data_test::perform() throw(libtest::test_exception) {

    test_1a();
    test_1b();
}


using namespace arma;


void orbital_data_test::test_1a() throw(libtest::test_exception) {

    static const char *testname = "orbital_data_test::test_1a()";

    // Perform operation
    try {

    // Test diagonalization of a restricted and symmetric density matrix

    size_t nao = test00_data::k_nao;
    size_t nmo = test00_data::k_nmo;
    test00_data data;

    mat s(nao, nao), c(nao, nmo);
    read_matrix(data, testname, "s", s);
    read_matrix(data, testname, "c_a", c);

    mat dm = randu(nmo, nmo);
    dm = c * (dm + dm.t()) * c.t();

    orbital_data orb(s, c, dm);

    // Check result
    const mat &u = orb.get_coeff();
    const vec &ev = orb.get_occ();
    if (accu(abs(u * diagmat(ev) * u.t() - dm) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform (1).");
    }
    if (accu(abs(u.t() * s * dm * s * u - diagmat(ev)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform (2).");
    }
    if (accu(abs(u.t() * s * u - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform (3).");
    }

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


void orbital_data_test::test_1b() throw(libtest::test_exception) {

    static const char *testname = "orbital_data_test::test_1b()";

    // Test diagonalization of a restricted and symmetric density matrix

    try {

    size_t nao = test00_data::k_nao;
    size_t nmo = test00_data::k_nmo;
    test00_data data;

    mat s(nao, nao), c(nao, nmo);
    read_matrix(data, testname, "s", s);
    read_matrix(data, testname, "c_b", c);

    mat dm = randu(nmo, nmo);
    dm = c * (dm + dm.t()) * c.t();

    // Perform operation

    orbital_data orb(s, c, dm);

    // Check result
    const mat &u = orb.get_coeff();
    const vec &ev = orb.get_occ();
    if (accu(abs(u * diagmat(ev) * u.t() - dm) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform (1).");
    }
    if (accu(abs(u.t() * s * dm * s * u - diagmat(ev)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform (2).");
    }
    if (accu(abs(u.t() * s * u - eye(nmo, nmo)) > 1e-12) != 0) {
        fail_test(testname, __FILE__, __LINE__, "Bad transform (3).");
    }

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


} // namespace libwfa
