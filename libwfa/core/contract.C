#include "contract.h"

namespace libwfa{
using namespace arma;

double contract::perform(const ab_matrix &tdm, const ab_matrix &om,
        const std::string op1, const std::string op2, const char spin) const{

/** Test to see if the input is possible
if(!((((op1=="x")#(op1=="y")#(op1=="z")#(op1=="s"))&&((op2=="x")#(op2=="y")#
    (op2=="z"))#(op2=="s"))#(((op1=="xx")#(op1=="yy")#(op1=="zz")#(op1=="s"))&&
    ((op2=="xx")#(op2=="yy")#(op2=="zz")#(op2=="s")))))
        {
        return 0;
        }//Endif
**/

    switch (spin){
        case 'a':{
                double omval_a=accu(om.alpha());
                double expval_a=0;
                Mat<double> exp_a;
                Mat<double> m1;
                Mat<double> m2;

                if (op1=="x") m1=m_x;
                if (op1=="y") m1=m_y;
                if (op1=="z") m1=m_z;
                if (op1=="xx") m1=m_xx;
                if (op1=="yy") m1=m_yy;
                if (op1=="zz") m1=m_zz;
                if (op1=="s") m1=m_s;

                if (op2=="x") m2=m_x;
                if (op2=="y") m2=m_y;
                if (op2=="z") m2=m_z;
                if (op2=="xx") m2=m_xx;
                if (op2=="yy") m2=m_yy;
                if (op2=="zz") m2=m_zz;
                if (op2=="s") m2=m_s;

                exp_a=(tdm.alpha()*m1)%(m2*tdm.alpha());
                expval_a=accu(exp_a)/omval_a;
                return expval_a;
                break;
         }//endcase

        case 'b':{
            double omval_b=accu(om.beta());
            double expval_b=0;
            Mat<double> exp_b;
            Mat<double> m1;
            Mat<double> m2;

            if (op1=="x") m1=m_x;
            if (op1=="y") m1=m_y;
            if (op1=="z") m1=m_z;
            if (op1=="xx") m1=m_xx;
            if (op1=="yy") m1=m_yy;
            if (op1=="zz") m1=m_zz;
            if (op1=="s") m1=m_s;

            if (op2=="x") m1=m_x;
            if (op2=="y") m1=m_y;
            if (op2=="z") m1=m_z;
            if (op2=="xx") m1=m_xx;
            if (op2=="yy") m1=m_yy;
            if (op2=="zz") m1=m_zz;
            if (op2=="s") m1=m_s;

            exp_b=(tdm.beta()*m1)%(m2*tdm.beta());
            expval_b=accu(exp_b)/omval_b;
            return expval_b;
            break;
        }//endcase
            default:
                return 0;

        } //end switch

}//end fct











}//end namespace libwfa
