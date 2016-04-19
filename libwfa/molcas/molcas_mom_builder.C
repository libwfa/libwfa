#include <libwfa/libwfa_exception.h>
#include "molcas_mom_builder.h"

namespace libwfa {

using namespace arma;

double molcas_mom_builder::perform(const mat &dm,
    size_t c1, size_t n1, size_t c2, size_t n2) const{

    const mat &m1 = get(c1, n1);
    const mat &m2 = get(c2, n2);

    return accu((dm * m1) % (m2 * dm));
}


double molcas_mom_builder::perform (const mat &dm, size_t c, size_t n) const{

    const mat &m = get(c, n);

    return trace(dm * m);
}

mat &molcas_mom_builder::set(size_t c, size_t n) {

    size_t pos = get_index(c, n);
    if (pos > m_op.size()) {
        throw libwfa_exception("mom_builder", "get(size_t, size_t)",
                __FILE__, __LINE__, "n");
    }

    m_op[pos] = new mat(m_nao, m_nao);

    return *m_op[pos];
}

const mat &molcas_mom_builder::get(size_t c, size_t n) const{

    size_t pos = get_index(c, n);
    if (pos > m_op.size()) {
        throw libwfa_exception("mom_builder", "get(size_t, size_t)",
                __FILE__, __LINE__, "n");
    }

    if (m_op[pos] == 0)
    {
        std::cout << "Matrix not set, c = " << c << ", n = " << n;
        throw libwfa_exception("mom_builder", "get",
                __FILE__, __LINE__, "Matrix not set");
    }

    return *m_op[pos];
}

void molcas_mom_builder::initialize(size_t pos) const {
    
    if (m_op[pos] != 0) return;
    
    throw libwfa_exception("mom_builder", "init",
                __FILE__, __LINE__, "not implemented");
       // hier weiter
/*        H5::H5std_string key;
        switch (pos) {
            case 1: key = "AO_MLTPL_X";
            case 2: key = "AO_MLTPL_Y";
            case 3: key = "AO_MLTPL_Z";
            case 4: key = "AO_MLTPL_XX";
            case 5: key = "AO_MLTPL_YY";
            case 6: key = "AO_MLTPL_ZZ";
        }*/
     
}

}//end namespace libwfa