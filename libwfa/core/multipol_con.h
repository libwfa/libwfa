#ifndef LIBWFA_MULTIPOL_CON_H
#define LIBWFA_MULTIPOL_CON_H

#include "multipol_con_i.h"

namespace libwfa{

using namespace arma;


/** \brief Implementation of multipol_con_i

	\ingroup libwfa
 **/
class multipol_con : public multipol_con_i {
private:
    const arma::Mat<double> &m_x; //!< Matrix for multipol moment \f$r_x\f$
    const arma::Mat<double> &m_xx;//!< Matrix for multipol moment \f$r²_x\f$
    const arma::Mat<double> &m_y;//!< Matrix for multipol moment \f$r_y\f$
    const arma::Mat<double> &m_yy;//!< Matrix for multipol moment \f$r²_y\f$
    const arma::Mat<double> &m_z;//!< Matrix for multipol moment \f$r_z\f$
    const arma::Mat<double> &m_zz;//!< Matrix for multipol moment \f$r²_z\f$
    const arma::Mat<double> &m_s;//!< Overlap matrix


public:
    /**\brief Constructor, setting needed matrices
       \param mx Matrix for multipol moment \f$r_x\f$
       \param mxx Matrix for multipol moment \f$r²_x\f$
       \param my Matrix for multipol moment \f$r_y\f$
       \param myy Matrix for multipol moment \f$r²_y\f$
       \param mz Matrix for multipol moment \f$r_z\f$
       \param mzz Matrix for multipol moment \f$r²_z\f$
       \param ms Overlap matrix

     **/
    multipol_con(
		const arma::Mat<double> &mx, const arma::Mat<double> &mxx,
		const arma::Mat<double> &my, const arma::Mat<double> &myy,
		const arma::Mat<double> &mz, const arma::Mat<double> &mzz,
        const arma::Mat<double> &ms) :
		m_x(mx), m_xx(mxx), m_y(my), m_yy(my), m_z(mz), m_zz(mzz), m_s(ms) {

    }
    /** \brief Performs the computation
                 \f$ \mathrm{Tr}\!\left[\gamma O_1 \gamma O_2 \right] \f$
        \param dm Density matrix \f$\gamma\f$
        \param op1 First operator \f$ O_1\f$
        \param op2 Second operator \f$ O_2\f$
        \return result
     **/
    virtual double perform (const Mat<double> &dm, const std::string &op1,
    		const std::string &op2) const;
    /**\brief Performs the computation
                     \f$ \mathrm{Tr}\!\left[\gamma O\right] \f$
      \param dm Density matrix \f$\gamma\f$
      \param op Operator \f$ O\f$
      \return result
     **/

    virtual double perform (const Mat<double> &dm, const std::string &op) const;

private:
    /**\brief Gives back the desired operator as a matrix
            \param op Operator
     **/
    const Mat<double> &retrieve_op(const std::string &op) const;
};


}//end namespace libwfa

#endif // LIBWFA_MULTIPOL_CON_H
