#ifndef LIBWFA_POP_ANALYSIS_I_H
#define LIBWFA_POP_ANALYSIS_I_H

#include <armadillo>
#include <vector>

namespace libwfa {

/** \brief Interface for classes to perform a certain population analysis

    Simple interface for population analysis.

    \ingroup libwfa
 **/
class pop_analysis_i {
public:
    virtual ~pop_analysis_i() { }

    /** \brief Perform population analysis for given density
        \param[in] d_bb Density matrix in AO basis
        \param[out] p Resulting population analysis (length of p = # atoms)

        It is assumed that the function deletes the data in p.
     **/
    virtual void perform(
            const arma::Mat<double> &d_bb, std::vector<double> &p) = 0;
};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_I_H
