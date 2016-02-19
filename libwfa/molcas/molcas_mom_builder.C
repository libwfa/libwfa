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

const mat &molcas_mom_builder::get(size_t c, size_t n) const{

    size_t pos = get_index(c, n);
    if (pos > m_op.size()) {
        throw libwfa_exception("mom_builder", "get(size_t, size_t)",
                __FILE__, __LINE__, "n");
    }

    initialize(pos);

    return *m_op[pos];
}

void molcas_mom_builder::initialize(size_t pos) const {
    
    if (m_op[pos] != 0) return;
    
}

}//end namespace libwfa