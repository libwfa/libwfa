#ifndef LIBWFA_AB_ORBITAL_SELECTOR_H
#define LIBWFA_AB_ORBITAL_SELECTOR_H

#include "ab_object.h"
#include "orbital_selector.h"

namespace libwfa {

/** \brief Container for spin-selectors

    \ingroup libwfa
 **/
class ab_orbital_selector : public ab_object<orbital_selector> {
public:
    /** \brief Default constructor
        \param aeqb If true, alpha == beta spin vector
     **/
    ab_orbital_selector(bool aeqb = false) :
        ab_object<orbital_selector>(aeqb) { }


    /** \brief Constructor for alpha == beta
        \param nidx Number of indexes in range
     **/
    ab_orbital_selector(size_t nidx) : ab_object<orbital_selector>(true) {
        alpha() = orbital_selector(nidx);
    }

    /** \brief Constructor for alpha != beta
        \param nidx_a Number of alpha-spin indexes
        \param nidx_b Number of beta-spin indexes
     **/
    ab_orbital_selector(size_t nidx_a, size_t nidx_b) :
        ab_object<orbital_selector>(false) {
        alpha() = orbital_selector(nidx_a);
        beta()  = orbital_selector(nidx_b);
    }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nidx_a() const { return alpha().n_indexes(); }

    /** \brief Return the number of beta-spin rows
     **/
    size_t nidx_b() const { return beta().n_indexes(); }
};

} // namespace libwfa

#endif // LIBWFA_AB_ORBITAL_SELECTOR_H
