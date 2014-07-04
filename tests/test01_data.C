#include "test01_data.h"

namespace libwfa {

using namespace arma;

const size_t test01_data::k_nao = 6;
const size_t test01_data::k_nmo = 6;
const size_t test01_data::k_natoms = 2;
const size_t test01_data::k_nstat = 2;

test01_data::test01_data() : test_data_base("test01"),
    m_atnum(k_natoms), m_nch(k_natoms), m_bf2nuc(k_nao, fill::ones),
    m_popref(k_natoms, k_nstat+1, fill::zeros) {

    m_atnum(0) = 2;
    m_atnum(1) = 3;
    m_nch(0) = 2.0;
    m_nch(1) = 3.0;
    m_bf2nuc(0) = 0;
}


} // namespace libwfa
