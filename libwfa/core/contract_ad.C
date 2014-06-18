#include "contract_ad.h"

namespace libwfa{
using namespace arma;


double contract_ad::perform (const ab_matrix &ad,
            const double prom, const std::string op1, const char spin) const{
    switch (spin){
        case 'a':{
                double expval_a=0;
                Mat<double> exp_a;
                Mat<double> m1;

                if (op1=="x") m1=m_x;
                if (op1=="y") m1=m_y;
                if (op1=="z") m1=m_z;
                if (op1=="xx") m1=m_xx;
                if (op1=="yy") m1=m_yy;
                if (op1=="zz") m1=m_zz;
                if (op1=="s") m1=m_s;

                exp_a=ad.alpha()*m1;
                expval_a=trace(exp_a)/prom;
                return expval_a;
                break;
         }

        case 'b':{
                double expval_b=0;
                Mat<double> exp_b;
                Mat<double> m1;

                if (op1=="x") m1=m_x;
                if (op1=="y") m1=m_y;
                if (op1=="z") m1=m_z;
                if (op1=="xx") m1=m_xx;
                if (op1=="yy") m1=m_yy;
                if (op1=="zz") m1=m_zz;
                if (op1=="s") m1=m_s;

                exp_b=ad.beta()*m1;
                expval_b=trace(exp_b)/prom;
                return expval_b;
                break;
        }
            default:
                return 0;

        } //end switch
}//end fct

}//end namespace libwfa
