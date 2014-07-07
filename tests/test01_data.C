#include "test01_data.h"

namespace libwfa {

using namespace arma;

const size_t test01_data::k_nao = 6;
const size_t test01_data::k_nmo = 6;
const size_t test01_data::k_natoms = 2;
const size_t test01_data::k_nstat = 2;

test01_data::test01_data() : test_data_base("test01"),
    m_atnum(k_natoms), m_nch(k_natoms), m_bf2nuc(k_nao, fill::ones),
    m_popref_a(k_natoms, k_nstat+1), m_popref_b(k_natoms, k_nstat+1) {

    m_atnum(0) = 2;
    m_atnum(1) = 3;
    m_nch(0) = 2.0;
    m_nch(1) = 3.0;
    m_bf2nuc(0) = 0;
    
    m_popref_a(0,0) = -0.999587924698302;
    m_popref_a(1,0) = -2.000412075301708;

    m_popref_a(0,1) = -0.999420401906908;
    m_popref_a(1,1) = -2.0005795980931;

    m_popref_a(0,2) = -0.999420401906908;
    m_popref_a(1,2) = -2.000579598093101;

    // adjust data
    m_popref_b(0,0) = -0.999422369360491;  
    m_popref_b(1,0) = -1.000577630639517;

    m_popref_b(0,1) = -0.999419229022855;
    m_popref_b(1,1) = -1.000580770977153;

    m_popref_b(0,2) = -0.999419229022855;
    m_popref_b(1,2) = -1.000580770977153;
               
}


} // namespace libwfa
