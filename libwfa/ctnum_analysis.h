#ifndef LIBWFA_CTNUM_ANALYSIS_H
#define LIBWFA_CTNUM_ANALYSIS_H

#include <vector>
#include "ctnum_analysis_i.h"

namespace libwfa {

/** \brief Charge transfer (CT) number analysis.

    Main implementation of CT number analysis.

    TODO: rename if there are more variants of CT number analysis to come...

    \ingroup libwfa
 **/
class ctnum_analysis : public ctnum_analysis_i {
private:
    size_t m_natoms; //!< Number of atoms
    const std::vector<size_t> &m_b2c; //!< Map of basis functions to atoms

public:
    /** \brief Constructor
        \param b2c Map of atomic basis functions to atoms/nuclei
     **/
    ctnum_analysis(const std::vector<size_t> &b2c);

    virtual ~ctnum_analysis() { }

    /** \copydoc ctnum_analysis_i::dim
     **/
    virtual size_t dim() const {
        return m_natoms;
    }

    /** \copydoc ctnum_analysis_i::perform
     **/
    virtual void perform(const arma::Mat<double>& om_ao,
            arma::Mat<double> &om_at);
};

} // namespace libwfa

#endif
