#include "test02_data.h"

namespace libwfa {

using namespace arma;

const size_t test02_data::k_nao = 13;
const size_t test02_data::k_nmo = 13;
const size_t test02_data::k_natoms = 3;
const size_t test02_data::k_nstat = 2;

test02_data::test02_data() : test_data_base("test02"),
    m_atnum(k_natoms), m_nch(k_natoms), m_bf2nuc(k_nao, fill::zeros),
    m_popref(k_natoms, k_nstat+1) {

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
    
    m_popref(0,0) = -3.346223094764419;
    m_popref(1,0) = -0.326888452617794;
    m_popref(2,0) = -0.326888452617789;
    
    m_popref(0,1) = -2.984590364168055;
    m_popref(1,1) = -0.507704817915983;
    m_popref(2,1) = -0.507704817915964;
    
    m_popref(0,2) = -3.000874165011687;
    m_popref(1,2) = -0.499562917494166;
    m_popref(2,2) = -0.499562917494148;
}


} // namespace libwfa
