#include <sstream>
#include "compare_ref.h"

namespace libwfa {

using namespace arma;


void compare_ref::compare(const char *test, Mat<double> &m, Mat<double> &m_ref,
    double thresh) throw(libtest::test_exception) {

    if (m.n_cols != m_ref.n_cols) {
        throw libtest::test_exception("compare_ref::compare()",
            __FILE__, __LINE__, "# cols");
    }
    if (m.n_rows != m_ref.n_rows) {
        throw libtest::test_exception("compare_ref::compare()",
            __FILE__, __LINE__, "# rows");
    }

    for (size_t i = 0; i < m.n_rows; i++) {
        for (size_t j = 0; j < m.n_cols; j++) {
            if (fabs(m(i, j) - m_ref(i, j)) > thresh) {
                std::ostringstream ess;
                ess << "In " << test << ": ";
                ess << "Result does not match reference at element ("
                        << i << ", " << j << "): " << m(i, j)
                        << " (act) vs " << m_ref(i, j) << " (ref), "
                        << m(i, j) - m_ref(i, j) << " (diff)";
                throw libtest::test_exception("compare_ref::compare()",
                        __FILE__, __LINE__, ess.str().c_str());
            }
        }
    }
}


void compare_ref::compare(const char *test, ab_matrix &m, ab_matrix &m_ref,
    double thresh) throw(libtest::test_exception) {

    if (m.is_alpha_eq_beta() != m_ref.is_alpha_eq_beta()) {
        throw libtest::test_exception("compare_ref::compare()",
            __FILE__, __LINE__, "a != b");
    }
    compare_ref::compare(test, m.alpha(), m_ref.alpha(), thresh);
    compare_ref::compare(test, m.beta(), m_ref.beta(), thresh);
}


} // namespace libwfa
