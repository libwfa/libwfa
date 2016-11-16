#ifndef LIBWFA_DENS_MOM_H
#define LIBWFA_DENS_MOM_H

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/mom_builder_i.h>

namespace libwfa {

using namespace arma;

/** \brief Simple multipole analysis of the wavefunctions.

    Computes the usual dipole and quadrupole moments.
    The information is combined to compute the size of the
    electron distribution.

    TODO: Analysis of spin, quadrupole moments, ellipticity, ...

    \ingroup libwfa
 **/
class dens_mom {
private:
    size_t m_nmax; //!< Highest moment
    double m_nelec;   //!< Trace of the DMAT
    arma::vec m_dip;  //!< Electronic dipole moment
    arma::vec m_quad; //!< Electronic quadrupole moment (only diagonal elements for now)
    double m_tnuc;       //!< Total nuclear charge
    arma::vec m_ndip;  //!< Nuclear dipole moment
    arma::vec m_nquad; //!< Nuclear quadruople moment (diagonal elements)

public:
    /** \brief Constructor
        \param bld Moment builder
        \param sdm Density matrix
        \param maxmm Highest moment to perform analysis
     **/
    dens_mom(const mom_builder_i &bld,
        const ab_matrix &sdm, const arma::mat &coordinates,
        const arma::vec &atomic_charges, size_t maxmm = 2);

    /** \brief Destructor
     **/
    virtual ~dens_mom() { }

    void analyse(std::ostream &out, size_t offset = 2) const;

private:
    virtual void print_header(std::ostream &out, size_t off) const;

    virtual void analysis(std::ostream &out, size_t off) const;

    void calculate(const mom_builder_i &bld, const ab_matrix &sdm,
        const arma::mat &coordinates, const arma::vec &atomic_charges);

    /** \brief Print armadillo vector
        \param out Output stream
        \param vec Vector to print
        \param width Width of each element
     **/
    static void print(std::ostream &out,
        const arma::vec &vec, size_t width = 10);};

} // namespace libwfa


#endif
