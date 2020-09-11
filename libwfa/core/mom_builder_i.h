#ifndef LIBWFA_MOM_BUILDER_I_H
#define LIBWFA_MOM_BUILDER_I_H

#include <armadillo>
#include <string>

namespace libwfa {

/** \brief Interface for builders of spatial moments of a density matrix

    \ingroup libwfa
 **/
class mom_builder_i {
public:
    virtual ~mom_builder_i() {}

    /** \brief Return the maximum moment that can be constructed
     **/
    virtual size_t max_moment() const = 0;

    /** \brief Performs the computation
            \f$ \mathrm{Tr}\!\left[\gamma O_1 \gamma O_2 \right] \f$
        \param dm Density matrix \f$\gamma\f$
        \param c1 Direction of first operator (0 = x, 1 = y, 2 = z)
        \param n1 Exponent of first operator
        \param c2 Direction of second operator (0 = x, 1 = y, 2 = z)
        \param n2 Exponent of second operator
        \return Resulting moment

        Computes the moment of the transition density matrix according to the
        formula above. Each operator \f$ O_i \f$ is determined by \f$ c_i \f$
        and \f$ n_i \f$ where \f$ n_i \f$ is the exponent and \f$ c_i \f$ is
        the cartesian direction, e.g. \f$ (c_i, n_i) = (0,2) \f$ refers to
        the operator \f$ x^2 \f$ .

        TODO: extend to mixed direction moments
     **/
    virtual double perform(const arma::mat &dm,
        size_t c1, size_t n1, size_t c2, size_t n2) const = 0;

    /** \brief Performs the computation
            \f$ \mathrm{Tr}\!\left[\gamma O\right] \f$

        \param dm Density matrix \f$\gamma\f$
        \param c Direction of operator (0 = x, 1 = y, 2 = z)
        \param n Exponent of operator
        \return Resulting moment

        Computes the moment of the density matrix according to the
        formula above. The operator \f$ O \f$ is determined by \f$ c \f$
        and \f$ n \f$ where \f$ n \f$ is the exponent and \f$ c \f$ is
        the cartesian direction, e.g. \f$ (c, n) = (1,1) \f$ refers to the
        operator \f$ y \f$ .
     **/
    virtual double perform(const arma::mat &dm, size_t c, size_t n) const = 0;
};

} // namespace libwfa


#endif // LIBWFA_MULTIPOL_CON_I_H
