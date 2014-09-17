#ifndef LIBWFA_EXCITON_MOMENTS_H
#define LIBWFA_EXCITON_MOMENTS_H

#include <cstddef>
#include <armadillo>
#include <libwfa/core/ab_object.h>
#include <libwfa/libwfa_exception.h>


namespace libwfa {

using namespace arma;

/** \brief Container of exciton moments

    Stores the moments of an electron-hole pair

    \ingroup libwfa
 **/
class exciton_moments {
private:
    size_t m_nmax; //!< Highest moment
    arma::mat m_mom; //!< Computed moments

public:
    /** \brief Default constructor
     **/
    exciton_moments(size_t nmax = 1);

    /** \brief Set data for specific moment
     **/
    void set(size_t ne, size_t nh, const arma::vec &m) {
        m_mom.col(determine_loc(ne, nh)) = m;
    }

    arma::subview_col<double> get(size_t ne, size_t nh) const {
        return m_mom.col(determine_loc(ne, nh));
    }

private:
    size_t determine_loc(size_t ne, size_t nh) const;

}; // class exciton_moments


/** \brief Exciton moments with spin

    \ingroup libwfa
 **/
class ab_exciton_moments : public ab_object<exciton_moments> {
public:
    /** \brief Default constructor
        \param nmax Maximum moment (default: 2)
        \param aeqb \f$ \alpha = \beta \f$
     **/
    ab_exciton_moments(size_t nmax = 2, bool aeqb = false);

}; // class ab_exciton_moments


inline exciton_moments::exciton_moments(size_t nmax) : m_nmax(nmax),
    m_mom(3, ((nmax + 1) * (nmax + 2)) / 2 - 1, arma::fill::zeros) {

}


inline size_t exciton_moments::determine_loc(size_t ne, size_t nh) const {

    size_t n = ne + nh;
#ifdef LIBWFA_DEBUG
    if (n > m_nmax) {
        throw libwfa_exception("exciton_moments",
            "determine_loc(size_t, size_t)", __FILE__, __LINE__, "N/A");
    }
#endif

    return (n * (n + 1)) / 2 + ne - 1;
}


inline ab_exciton_moments::ab_exciton_moments(size_t nmax, bool aeqb) :
    ab_object<exciton_moments>(aeqb) {

    alpha() = exciton_moments(nmax);
    if (! aeqb) {
        beta() = exciton_moments(nmax);
    }
}



} // namespace libwfa


#endif
