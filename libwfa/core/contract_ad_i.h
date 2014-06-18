#ifndef CONTRACT_AD_I_H_
#define CONTRACT_AD_I_H_
#include <string>

namespace libwfa{

class ab_matrix;

class contract_ad_i{
public:

    virtual ~contract_ad_i() {}

    virtual double perform(const ab_matrix &ad,
            const double prom, const std::string op1, const char spin) const=0;

};

}//end namespace libwfa



#endif /* CONTRACT_AD_I_H_ */
