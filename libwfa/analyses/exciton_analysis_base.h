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
        \param offset Line offset
     **/
    void analyse(std::ostream &out, size_t offset = 2) const;

protected:
    /** \brief Write analysis header to output stream
        @param out Output stream
        @param offset Line offset
     **/
    virtual void print_header(std::ostream &out, size_t offset) const = 0;

    /** \brief Analyse exciton moments and print results
        @param out Output stream
        @param mom Exciton moments
        @param offset Line offset
     **/
    virtual void analysis(std::ostream &out,
            const exciton_moments &mom, size_t offset) const = 0;

    /** \brief Print armadillo vector
        \param out Output stream
        \param vec Vector to print
        \param width Width of each element
     **/
    static void print(std::ostream &out,
        const arma::vec &vec, size_t width = 10);

    /** \brief Form total exciton moments
        \param a Alpha exciton moment
        \param b Beta exciton moment
        \param res Total exciton moment
     **/
    virtual void combine(const exciton_moments &a, const exciton_moments &b,
        exciton_moments &res) const = 0;
};

} // namespace libwfa


#endif
