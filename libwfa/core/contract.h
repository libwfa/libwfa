#ifndef CONTRACT_H_
#define CONTRACT_H_

#include "ab_matrix.h"
#include "contract_i.h"

namespace libwfa{

using namespace arma;

class contract : public contract_i {
private:
    const Mat<double> &m_x;
    const Mat<double> &m_xx;
    const Mat<double> &m_y;
    const Mat<double> &m_yy;
    const Mat<double> &m_z;
    const Mat<double> &m_zz;
    const Mat<double> &m_s;


public:

    contract (const Mat<double> &mx, const Mat<double> &mxx, const Mat<double> &my,
        const Mat<double> &myy, const Mat<double> &mz, const Mat<double> &mzz,
        const Mat<double> &ms):m_x(mx), m_xx(mxx), m_y(my), m_yy(my), m_z(mz),
        m_zz(mzz), m_s(ms){}

    virtual double perform (const ab_matrix &tdm, const ab_matrix &om,
            const std::string op1,const std::string op2,
            const char spin) const;//endfct
};


}//end namespace libwfa

#endif /* CONTRACT_H_ */
