#ifndef LIBWFA_POP_LOEWDIN_H
#define LIBWFA_POP_LOEWDIN_H

#include <cstddef>
#include "pop_analysis_i.h"

namespace libwfa {

/** \brief Implementation of pop_analysis_i to perform LOEWDIN population
        analysis

    \ingroup libwfa
**/
class pop_loewdin : public pop_analysis_i {
private:
    size_t m_nparts; //!< Number of atoms
    arma::Mat<double> m_sh; //!< Square root of overlap matrix
    const arma::Col<size_t> &m_b2p; //!< Map of basis functions to parts

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param b2p Map of atomic basis functions to molecular parts
        \param p0 Population data to add
     **/
    pop_loewdin(const arma::Mat<double> &s, const arma::Col<size_t> &b2p);

    /** \brief Destructor
     **/
    virtual ~pop_loewdin() { }

    /** \copydoc pop_analysis_i::size
     **/
    virtual size_t size() const {  return m_nparts; }

    /** \copydoc pop_analysis_i::perform
     **/
    virtual void perform(const arma::Mat<double> &d_bb,
            arma::Col<double> &p) const;

};


} // namespace libwfa

#endif // LIBWFA_POP_LOEWDIN_H
