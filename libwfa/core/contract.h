#ifndef LIBWFA_CONTRACT_H
#define LIBWFA_CONTRACT_H

#include "contract_i.h"

namespace libwfa{

using namespace arma;


/** \brief Implementation of contract_i

  	TODO: Rename (more concrete naming to avoid confusion with general contractions)
  	TODO: Document properly (in particular the constructor).

	\ingroup libwfa
 **/
class contract : public contract_i {
private:
    const arma::Mat<double> &m_x;
    const arma::Mat<double> &m_xx;
    const arma::Mat<double> &m_y;
    const arma::Mat<double> &m_yy;
    const arma::Mat<double> &m_z;
    const arma::Mat<double> &m_zz;
    const arma::Mat<double> &m_s;


public:
    contract(
		const arma::Mat<double> &mx, const arma::Mat<double> &mxx,
		const arma::Mat<double> &my, const arma::Mat<double> &myy,
		const arma::Mat<double> &mz, const arma::Mat<double> &mzz,
        const arma::Mat<double> &ms) :
		m_x(mx), m_xx(mxx), m_y(my), m_yy(my), m_z(mz), m_zz(mzz), m_s(ms) {

    }

    virtual double perform (const Mat<double> &dm, const std::string &op1,
    		const std::string &op2) const;

    virtual double perform (const Mat<double> &dm, const std::string &op) const;

private:
    const Mat<double> &retrieve_op(const std::string &op) const;
};


}//end namespace libwfa

#endif // LIBWFA_CONTRACT_H
