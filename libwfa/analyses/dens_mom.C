#include <iomanip>
#include <libwfa/core/constants.h>
#include "dens_mom.h"

namespace libwfa {

using namespace arma;

dens_mom::dens_mom(const mom_builder_i &bld,
        const ab_matrix &sdm, const arma::mat &coordinates,
        const arma::vec &atomic_charges, size_t maxmm) :
            m_nmax(maxmm), m_dip(3), m_quad(3), m_ndip(3), m_nquad(3) {

    calculate(bld, sdm, coordinates, atomic_charges);
}

void dens_mom::analyse(std::ostream &out, size_t offset) const {
    print_header(out, offset);

    analysis(out, offset + 2);
}

void dens_mom::print_header(std::ostream &out, size_t off) const {

    out << std::string(off, ' ');
    out << "Multipole moment analysis of the density matrix" << std::endl;
}

void dens_mom::analysis(std::ostream &out, size_t off) const {
    std::string os(off, ' ');

    { // Dipole moment
        out << os << "Molecular charge:     " << std::string(17, ' ') << m_tnuc - m_nelec << std::endl;
        out << os << "Number of electrons:  " << std::string(17, ' ') << m_nelec << std::endl;
        arma::vec pvec = m_dip * constants::au2ang / m_nelec;
        out << os << "  Center of electronic charge [Ang]: ";
        print(out, pvec);

        out << os << "Total nuclear charge: " << std::string(17, ' ') << m_tnuc << std::endl;
        pvec = m_ndip * constants::au2ang / m_nelec;
        out << os << "  Center of nuclear charge [Ang]:    ";
        print(out, pvec);

        pvec = (m_ndip - m_dip) * constants::au2D;
        out << os << "Dipole moment [D]: " << std::string(20, ' ') << norm(pvec) << std::endl;
        out << os << "  Cartesian components [D]: " << std::string(9, ' ');
        print(out, pvec);
    }

    /* TODO: adjust for nuclear quadrupole moment and compute full tensor

    out << os << "Quadrupole moment [a.u.]: ";
    print(out, m_quad);
    vec trl_quad = (3*m_quad - accu(m_quad)) / 2.;
    out << os << "Traceless quadrupole moment [a.u.]: ";
    print(out, trl_quad); */

    { // Molecular size
        vec d2 = m_quad / m_nelec - m_dip % m_dip / (m_nelec * m_nelec);
        d2  *= constants::au2ang * constants::au2ang;
        double d  = sqrt(accu(d2));
        out << os << "RMS size of the density [Ang]:" << std::string(7, ' ')
                << std::setw(10) << d << std::endl;
        out << os << "  Cartesian components [Ang]:" << std::string(8, ' ');
        print(out, sqrt(d2));
    }
    out << std::endl;
}

void dens_mom::calculate(const mom_builder_i &bld, const ab_matrix &sdm,
        const arma::mat &coordinates, const arma::vec &atomic_charges) {

    arma::mat sdm_tr = sdm.alpha() + sdm.beta();

    m_nelec = bld.perform(sdm_tr, 'x', 0);

    for (size_t j = 0; j < 3; j++) {
        m_dip(j)  = bld.perform(sdm_tr, j, 1);
        m_quad(j) = bld.perform(sdm_tr, j, 2);
    }

    m_tnuc = accu(atomic_charges);
    m_ndip = coordinates * atomic_charges;
}

void dens_mom::print(std::ostream &out,
        const vec &vec, size_t width) {

    out << "[";
    if (vec.n_rows > 0) out << std::setw(width) << vec(0);
    for (size_t i = 1; i < vec.n_rows; i++) {
        out << ", " << std::setw(width) << vec(i);
    }
    out << "]" << std::endl;
}

}// end namespace libwfa

