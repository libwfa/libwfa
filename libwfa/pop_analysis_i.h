#ifndef LIBWFA_POP_ANALYSIS_I_H
#define LIBWFA_POP_ANALYSIS_I_H

#include <armadillo>
#include <vector>

namespace libwfa {

/** \brief Interface for classes to perform a specific type of population
        analysis

    Simple interface for specific types of population analysis.

    \ingroup libwfa
 **/
class pop_analysis_i {
public:
    virtual ~pop_analysis_i() { }

    /** \brief Return the length of a population data set generated (e.g. # atoms)
     **/
    virtual size_t size() const = 0;

    /** \brief Perform population analysis for given density
        \param[in] dm Density matrix in AO basis
        \param[out] p Resulting population analysis

        Routine to perform the population analysis for one density matrix. The
        output vector should be adjusted to the correct size by the routine.
     **/
    virtual void perform(const arma::Mat<double> &dm,
            std::vector<double> &p) = 0;

};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_I_H
