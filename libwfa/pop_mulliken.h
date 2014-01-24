#ifndef LIBWFA_POP_MULLIKEN_H
#define LIBWFA_POP_MULLIKEN_H

#include <cstddef>
#include "pop_analysis_i.h"

namespace libwfa {

/** \brief Implementation of pop_analysis_i to perform Mulliken population
        analysis

    \ingroup libwfa
**/
class pop_mulliken : public pop_analysis_i {
private:
    size_t m_natoms; //!< Number of atoms
	const std::vector<size_t> &m_b2c; //!< Map of basis functions to atoms
	const arma::Mat<double> &m_s; //!< Overlap matrix

public:
    /** \brief Constructor
        \param b2c Map of atomic basis functions to atoms/nuclei
        \param s Overlap matrix
     **/
    pop_mulliken(const std::vector<size_t> &b2c, const arma::Mat<double> &s);

    virtual ~pop_mulliken() { }

    /** \copydoc pop_analysis_i::size
     **/
    virtual size_t size() const {
        return m_natoms;
    }

    /** \copydoc pop_analysis_i::perform
     **/
    virtual void perform(const arma::Mat<double> &d_bb,
            std::vector<double> &p) const;

};


} // namespace libwfa

#endif // LIBWFA_POP_MULLIKEN_H
