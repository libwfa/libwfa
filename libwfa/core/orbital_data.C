#include "orbital_data.h"

namespace libwfa {

using namespace arma;


orbital_data::orbital_data(const mat &s, const mat &c, const mat &dm) {

    // Orthogonalize density matrix via C' S DM S C
    mat cinv = s * c;
    m_coeff = cinv.t() * dm * cinv;

    // Diagonalize orthogonal version of density matrix
    mat evec;
    eig_sym(m_ev, evec, m_coeff);

    // Back-transform eigenvectors to non-orthogonal basis
    m_coeff = c * evec;
}


} // namespace libwfa

