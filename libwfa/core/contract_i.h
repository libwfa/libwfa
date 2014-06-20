#ifndef LIBWFA_CONTRACT_I_H
#define LIBWFA_CONTRACT_I_H

#include <armadillo>
#include <string>

namespace libwfa {

/** \brief Interface to evaluate multipole moments w.r.t a density matrix

  	\ingroup libwfa
 **/
class contract_i {
public:
    virtual ~contract_i() {}

    /** \brief Performs the computation
    		\f$ \mathrm{Tr}\!\left[\gamma O_1 \gamma O_2 \right] \f$

        \param dm Density matrix \f$\gamma\f$
        \param op1 First operator \f$ O_1\f$
        \param op2 Second operator \f$ O_2 \f$
        \return Result
     **/
    virtual double perform(const arma::Mat<double> &dm,
    		const std::string &op1, const std::string &op2) const = 0;

    /** \brief Performs the computation
    		\f$ \mathrm{Tr}\!\left[\gamma O\right] \f$

        \param dm Density matrix \f$\gamma\f$
        \param op Operator \f$ O\f$
        \return Result
     **/
    virtual double perform(const arma::Mat<double> &dm,
    		const std::string &op) const = 0;
};

}//end namespace libwfa


#endif // LIBWFA_CONTRACT_I_H
