#ifndef LIBWFA_MOM_BUILDER_H
#define LIBWFA_MOM_BUILDER_H

#include "mom_builder_i.h"

namespace libwfa{

using namespace arma;


/** \brief Simple implementation of mom_builder_i

    \ingroup libwfa
 **/
class mom_builder : public mom_builder_i {
private:
    const arma::mat &m_s;  //!< Overlap matrix (zero-th moment)
    const arma::mat &m_x;  //!< Moment matrix \f$ x \f$
    const arma::mat &m_y;  //!< Moment matrix \f$ y \f$
    const arma::mat &m_z;  //!< Moment matrix \f$ z \f$
    const arma::mat &m_xx; //!< Moment matrix \f$ x^2 \f$
    const arma::mat &m_yy; //!< Moment matrix \f$ y^2 \f$
    const arma::mat &m_zz; //!< Moment matrix \f$ z^2 \f$

public:

    /**\brief Constructor, setting needed matrices
       \param s  Overlap matrix (zero-th moment)
       \param x  Moment matrix \f$ x \f$
       \param y  Moment matrix \f$ y \f$
       \param z  Moment matrix \f$ z \f$
       \param xx Moment matrix \f$ x^2 \f$
       \param yy Moment matrix \f$ y^2 \f$
       \param zz Moment matrix \f$ z^2 \f$

     **/
    mom_builder(const arma::mat &s, const arma::mat &x,
        const arma::mat &y, const arma::mat &z,
        const arma::mat &xx, const arma::mat &yy,
        const arma::mat &zz) :
        m_s(s), m_x(x), m_y(y), m_z(z), m_xx(xx), m_yy(yy), m_zz(zz) { }

    /** \copydoc mom_builder_i::max_moment
     **/
    virtual size_t max_moment() const { return 2; }

    /** \copydoc mom_builder_i::perform
     **/
    virtual double perform(const arma::mat &dm,
        size_t c1, size_t n1, size_t c2, size_t n2) const;

    /** \copydoc mom_builder_i::perform
     **/
    virtual double perform(const arma::mat &dm,
        size_t c, size_t n) const;

private:
    /** \brief Returns the desired operator as a matrix
        \param c Coordinate direction
        \param n Exponent
        \return Operator matrix
     **/
    const mat &retrieve_op(size_t c, size_t n) const;
};


}//end namespace libwfa

#endif // LIBWFA_MOM_BUILDER_H
