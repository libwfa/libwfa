#ifndef LIBWFA_AB_VECTOR_H
#define LIBWFA_AB_VECTOR_H

#include <armadillo>
#include "ab_object.h"

namespace libwfa {

/** \brief Container for spin-vectors

    \ingroup libwfa
 **/
class ab_vector : public ab_object< arma::vec > {
public:
    /** \brief Default constructor
        \param aeqb If true, alpha == beta spin vector
     **/
    ab_vector(bool aeqb = false) : ab_object< arma::vec >(aeqb) { }

    /** \brief Constructor for alpha == beta
        \param nrows Number of rows
     **/
    ab_vector(size_t nrows) : ab_object< arma::vec >(true) {

        alpha() = arma::vec(nrows);
    }

    /** \brief Constructor for alpha != beta
        \param nrows_a Number of alpha-spin rows
        \param nrows_b Number of beta-spin rows
     **/
    ab_vector(size_t nrows_a, size_t nrows_b) :
        ab_object< arma::vec >(false) {

        alpha() = arma::vec(nrows_a);
        beta()  = arma::vec(nrows_b);
    }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nrows_a() const { return alpha().n_rows; }

    /** \brief Return the number of beta-spin rows
     **/
    size_t nrows_b() const { return beta().n_rows; }
};

} // namespace libwfa

#endif // LIBWFA_AB_VECTOR_H
