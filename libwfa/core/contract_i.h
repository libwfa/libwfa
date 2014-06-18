#ifndef CONTRACT_I_H_
#define CONTRACT_I_H_
#include <string>
namespace libwfa {

class ab_matrix;

class contract_i{
public:

    virtual ~contract_i() {}

    virtual double perform(const ab_matrix &tdm, const ab_matrix &om,
            const std::string op1, const std::string op2, const char spin) const=0;

};

}//end namespace libwfa


#endif /* CONTRACT_I_H_ */
