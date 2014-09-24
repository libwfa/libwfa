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
    size_t m_nparts; //!< Number of atoms
    const arma::mat &m_s; //!< Overlap matrix
    const arma::uvec &m_b2p; //!< Map of basis functions to parts

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param b2p Map of atomic basis functions to molecular parts
     **/
    pop_mulliken(const arma::mat &s, const arma::uvec &b2p);

    /** \brief Destructor
     **/
    virtual ~pop_mulliken() { }

    /** \copydoc pop_analysis_i::size
     **/
    virtual size_t size() const {  return m_nparts; }

    /** \copydoc pop_analysis_i::perform
     **/
    virtual void perform(const arma::mat &d_bb, arma::vec &p) const;

};


} // namespace libwfa

#endif // LIBWFA_POP_MULLIKEN_H
