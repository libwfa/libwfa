#include <armadillo>
#include <cmath>
#include <libwfa/export/cube_writer.h>
#include <libwfa/libwfa_exception.h>
#include "cube_writer_test.h"

using namespace arma;


namespace libwfa {

void cube_writer_test::perform() throw(libtest::test_exception) {

    test_1();
}


/** \brief Test construction of an empty grid
 **/
void cube_writer_test::test_1() {

    static const char *testname = "cube_writer_test::test_1()";

    try {

    unsigned int atnum[2] = { 8, 8 };
    double coord[6] = { 0.0, 0.0, -1.2, 0.0, 0.0, 1.2 };
    grid3d g(5, 0.5);
    cube_writer wr("cube_writer_test_1.cube", testname, "O2 test", g,
            Col<unsigned int>(atnum, 2), Mat<double>(coord, 3, 2));

    if (wr.npoints() != 125) fail_test(testname, __FILE__, __LINE__, "npts");

    Col<double> data(40);
    size_t i = 0, n = 0;
    for (; ! wr.complete(); i++) {
        data.randu();
        n += wr.write(data);
    }
    if (i != 4) fail_test(testname, __FILE__, __LINE__, "i");
    if (n != 125) fail_test(testname, __FILE__, __LINE__, "n");
    if (! wr.complete())
        fail_test(testname, __FILE__, __LINE__, "complete()");

    } catch (libwfa_exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }
}


} // namespace libwfa
