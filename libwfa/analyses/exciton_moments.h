#ifndef LIBWFA_EXCITON_MOMENTS_H
#define LIBWFA_EXCITON_MOMENTS_H

#include <cstddef>
#include <armadillo>
#include <libwfa/libwfa_exception.h>


namespace libwfa {

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

    /** \brief Highest moment
     **/
    size_t n_max() const { return m_nmax; }

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


} // namespace libwfa


#endif
