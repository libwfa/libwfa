#include "test02_data.h"

namespace libwfa {

using namespace arma;

const size_t test02_data::k_nao = 13;
const size_t test02_data::k_nmo = 13;
const size_t test02_data::k_natoms = 3;

test02_data::test02_data() : test_data_base("test02"),
    m_atnum(k_natoms), m_nch(k_natoms), m_bf2nuc(k_nao, fill::zeros) {

    m_atnum(0) = 8;
    m_atnum(1) = 1;
    m_atnum(2) = 1;
    m_nch(0) = 6.0;
    m_nch(1) = 1.0;
    m_nch(2) = 1.0;
    m_bf2nuc( 9) = 1;
    m_bf2nuc(10) = 1;
    m_bf2nuc(11) = 2;
    m_bf2nuc(12) = 2;
}


} // namespace libwfa
