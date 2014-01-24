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
    /** \brief Virtual destructor
     **/
    virtual ~pop_analysis_i() { }

    /** \brief Return the length of a population data generated (e.g. # atoms)
     **/
    virtual size_t size() const = 0;

    /** \brief Perform population analysis for given density
        \param[in] dm Density matrix in AO basis
        \param[out] p Resulting population data

        Routine to perform the population analysis for one density matrix. The
        routine should adjust the output vector to the correct size.
     **/
    virtual void perform(const arma::Mat<double> &dm,
            std::vector<double> &p) const = 0;

};

} // namespace libwfa

#endif // LIBWFA_POP_ANALYSIS_I_H
