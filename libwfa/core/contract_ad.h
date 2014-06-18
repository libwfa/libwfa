#ifndef CONTRACT_AD_H_
#define CONTRACT_AD_H_

#include "contract_ad_i.h"
#include "ab_matrix.h"

namespace libwfa{
using namespace arma;

class contract_ad : public contract_ad_i{
private:
    const Mat<double> &m_x;
    const Mat<double> &m_xx;
    const Mat<double> &m_y;
    const Mat<double> &m_yy;
    const Mat<double> &m_z;
    const Mat<double> &m_zz;
    const Mat<double> &m_s;

public:
    contract_ad (const Mat<double> &mx, const Mat<double> &mxx, const Mat<double> &my,
            const Mat<double> &myy, const Mat<double> &mz, const Mat<double> &mzz,
            const Mat<double> &ms):m_x(mx), m_xx(mxx), m_y(my), m_yy(my), m_z(mz),
            m_zz(mzz), m_s(ms){}

    virtual double perform (const ab_matrix &ad,
            const double prom, const std::string op1, const char spin) const;//endfct

};

}//end namespace libwfa


#endif /* CONTRACT_AD_H_ */
