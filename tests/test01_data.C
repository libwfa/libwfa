#include "test01_data.h"

namespace libwfa {

using namespace arma;

const size_t test01_data::k_nao = 6;
const size_t test01_data::k_nmo = 6;
const size_t test01_data::k_natoms = 2;
const size_t test01_data::k_nstat = 2;

test01_data::test01_data() : test_data_base("test01"),
    m_atnum(k_natoms), m_nch(k_natoms), m_bf2nuc(k_nao, fill::ones) {

    m_atnum(0) = 2;
    m_atnum(1) = 3;
    m_nch(0) = 2.0;
    m_nch(1) = 3.0;
    m_bf2nuc(0) = 0;

    double pma[6] = {
            -0.999587924698302, -2.000412075301708,
            -0.999420401906908, -2.0005795980931,
            -0.999420401906908, -2.000579598093101
    };
    double pmb[6] = {
            -0.999422369360491, -1.000577630639517,
            -0.999419229022855, -1.000580770977153,
            -0.999419229022855, -1.000580770977153
    };
    double pla[6] = {
            -0.999587924698302, -2.000412075301708,
            -0.999420401906908, -2.0005795980931,
            -0.999420401906908, -2.000579598093101
    };
    double plb[6] = {
            -0.999422369360491, -1.000577630639517,
            -0.999419229022855, -1.000580770977153,
            -0.999419229022855, -1.000580770977153
    };

    m_pop_mulliken[0] = Mat<double>(pma, k_natoms, k_nstat + 1);
    m_pop_mulliken[1] = Mat<double>(pmb, k_natoms, k_nstat + 1);
    m_pop_loewdin[0] = Mat<double>(pla, k_natoms, k_nstat + 1);
    m_pop_loewdin[1] = Mat<double>(plb, k_natoms, k_nstat + 1);
}


} // namespace libwfa
