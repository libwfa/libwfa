#ifndef LIBWFA_EXCITON_ANALYSIS_BASE_H
#define LIBWFA_EXCITON_ANALYSIS_BASE_H

#include <ostream>
#include "exciton_moments.h"


namespace libwfa {

/** \brief Base class for exciton analysis

    \ingroup libwfa
 **/
class exciton_analysis_base {
private:
    exciton_moments *m_mom[2]; //!< Computed exciton moments

public:
    /** \brief Constructor
        \param aeqb Is alpha-spin equal to beta-spin part
        \param maxmm Maximum moment (default: 2)
     **/
    exciton_analysis_base(bool aeqb, size_t maxmm = 2);

    /** \brief Destructor
     **/
    virtual ~exciton_analysis_base();

    /** \brief Is alpha equal to beta-spin
     **/
    bool is_alpha_eq_beta() const { return m_mom[1] == 0; }

    /** \brief Computed exciton moment
        \param spin If true: beta spin; else alpha spin
     **/
    exciton_moments &moment(bool spin) {
        return *m_mom[(spin && m_mom[1] ? 1 : 0)];
    }

    /** \brief Computed exciton moment
        \param spin If true: beta spin; else alpha spin
     **/
    const exciton_moments &moment(bool spin) const {
        return *m_mom[(spin && m_mom[1] ? 1 : 0)];
    }

    /** \brief Perform analysis
        \param out Output stream
     **/
    void analyse(std::ostream &out) const;

protected:
    virtual void print_header(std::ostream &out) const = 0;

    virtual void analysis(std::ostream &out,
            const exciton_moments &mom) const = 0;

    static void print(std::ostream &out,
        const arma::vec &vec, size_t width = 10);
};

} // namespace libwfa


#endif
