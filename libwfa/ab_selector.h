#ifndef LIBWFA_AB_SELECTOR_H
#define LIBWFA_AB_SELECTOR_H

#include "ab_object.h"
#include "selector.h"

namespace libwfa {

/** \brief Container for spin-selectors

    \ingroup libwfa
 **/
class ab_selector : public ab_object<selector> {
public:
    /** \brief Default constructor
        \param aeqb If true, alpha == beta spin vector
     **/
    ab_selector(bool aeqb = false) : ab_object<selector>(aeqb) { }


    /** \brief Constructor for alpha == beta
        \param nidx Number of indexes in range
     **/
    ab_selector(size_t nidx) : ab_object<selector>(true) {
        alpha() = selector(nidx);
    }

    /** \brief Constructor for alpha != beta
        \param nidx_a Number of alpha-spin indexes
        \param nidx_b Number of beta-spin indexes
     **/
    ab_selector(size_t nidx_a, size_t nidx_b) : ab_object<selector>(false) {
        alpha() = selector(nidx_a);
        beta()  = selector(nidx_b);
    }

    /** \brief Return the number of alpha-spin rows
     **/
    size_t nidx_a() const { return alpha().n_indexes(); }

    /** \brief Return the number of beta-spin rows
     **/
    size_t nidx_b() const { return beta().n_indexes(); }
};

} // namespace adcman

#endif // LIBWFA_AB_SELECTOR_H
