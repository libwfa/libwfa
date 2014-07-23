#include <libwfa/libwfa_exception.h>
#include "mom_builder.h"

namespace libwfa{


using namespace arma;

double mom_builder::perform(const Mat<double> &dm,
    size_t c1, size_t n1, size_t c2, size_t n2) const{

    const Mat<double> &m1 = retrieve_op(c1, n1);
    const Mat<double> &m2 = retrieve_op(c2, n2);

    return accu((dm * m1) % (m2 * dm));
}


double mom_builder::perform (const Mat<double> &dm, size_t c, size_t n) const{

    const Mat<double> &m = retrieve_op(c, n);

    return trace(dm * m);
}


const Mat<double> &mom_builder::retrieve_op(size_t c, size_t n) const{

    static const char clazz[] = "mom_builder";
    static const char method[] = "retrieve_op(size_t, size_t)";

    if (n == 0) return m_s;
    if (n > 2)
        throw libwfa_exception(clazz, method, __FILE__, __LINE__, "n");

    switch (c) {
    case 0: return (n == 1 ? m_x : m_xx);
    case 1: return (n == 1 ? m_y : m_yy);
    case 2: return (n == 1 ? m_z : m_zz);
    default:
        throw libwfa_exception(clazz, method, __FILE__, __LINE__, "c");
    }

}


}//end namespace libwfa
