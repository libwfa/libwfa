#include "test02_data.h"

namespace libwfa {

using namespace arma;

const size_t test02_data::k_nao = 13;
const size_t test02_data::k_nmo = 13;
const size_t test02_data::k_natoms = 3;
const size_t test02_data::k_nstat = 2;

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
    
    double pm[9] = {
            -3.346223094764419, -0.326888452617794, -0.326888452617789,
            -2.984590364168055, -0.507704817915983, -0.507704817915964,
            -3.000874165011687, -0.499562917494166, -0.499562917494148
    };
    double pl[9] = {
            -3.223819805200865, -0.388090097399570, -0.388090097399568,
            -2.956931836565553, -0.521534081717231, -0.521534081717219,
            -2.967265456772048, -0.516367271613983, -0.516367271613972
    };

    m_pop_mulliken = mat(pm, k_natoms, k_nstat + 1);
    m_pop_loewdin  = mat(pl, k_natoms, k_nstat + 1);
}


} // namespace libwfa
