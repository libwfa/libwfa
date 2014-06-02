#include "ex_analyse.h"

namespace libwfa {

using namespace arma;


double ex_analyse::ex_multip (const ab_matrix &tdm, const ab_matrix &op1,
        const ab_matrix &op2, const ab_matrix &om, char spin){

    switch (spin){
    case 'a':{
            const Mat<double> &tdm_a=tdm.alpha();
            const Mat<double> &op1_a=op1.alpha();
            const Mat<double> &op2_a=op2.alpha();
            double omval_a=accu(om.alpha());
            double expval_a=0;

            Mat<double> exp_a=(tdm_a*op2_a)%(op1_a*tdm_a);

            expval_a=accu(exp_a)/omval_a;
            return expval_a;
            break;
     }

    case 'b':{
            const Mat<double> &tdm_b=tdm.beta();
            const Mat<double> &op1_b=op1.beta();
            const Mat<double> &op2_b=op2.beta();
            double omval_b=accu(om.beta());
            double expval_b=0;

            Mat<double> exp_b=(tdm_b*op2_b)%(op1_b*tdm_b);

            expval_b=accu(exp_b)/omval_b;

            return expval_b;
            break;
    }
        default:
            return 0;

    } //end switch

}//end fct


void ex_analyse::ex_form (const ab_matrix &tdm, const ab_matrix &s,
                const ab_matrix &mxx, const ab_matrix &mx, const ab_matrix &myy,
                const ab_matrix &my, const ab_matrix &mzz,
                const ab_matrix &mz, const ab_matrix &om){

    rh[0][0]=ex_multip(tdm, mx, s, om, 'a');
    rh[1][0]=ex_multip(tdm, my, s, om, 'a');
    rh[2][0]=ex_multip(tdm, mz, s, om, 'a');

    re[0][0]=ex_multip(tdm, s, mx, om, 'a');
    re[1][0]=ex_multip(tdm, s, my, om, 'a');
    re[2][0]=ex_multip(tdm, s, mz, om, 'a');

    rh2[0][0]=ex_multip(tdm, mxx, s, om, 'a');
    rh2[1][0]=ex_multip(tdm, myy, s, om, 'a');
    rh2[2][0]=ex_multip(tdm, mzz, s, om, 'a');

    re2[0][0]=ex_multip(tdm, s, mxx, om, 'a');
    re2[1][0]=ex_multip(tdm, s, myy, om, 'a');
    re2[2][0]=ex_multip(tdm, s, mzz, om, 'a');

    rhre[0][0]=ex_multip(tdm, mx, mx, om, 'a');
    rhre[1][0]=ex_multip(tdm, my, my, om, 'a');
    rhre[2][0]=ex_multip(tdm, mz, mz, om, 'a');

    if (!tdm.is_alpha_eq_beta()){

        rh[0][1]=ex_multip(tdm, mx, s, om, 'b');
        rh[1][1]=ex_multip(tdm, my, s, om, 'b');
        rh[2][1]=ex_multip(tdm, mz, s, om, 'b');
        rh2[0][1]=ex_multip(tdm, mxx, s, om, 'b');
        rh2[1][1]=ex_multip(tdm, myy, s, om, 'b');
        rh2[2][1]=ex_multip(tdm, mzz, s, om, 'b');
        re[0][1]=ex_multip(tdm, s, mx, om, 'b');
        re[1][1]=ex_multip(tdm, s, my, om, 'b');
        re[2][1]=ex_multip(tdm, s, mz, om, 'b');
        re2[0][1]=ex_multip(tdm, s, mxx, om, 'b');
        re2[1][1]=ex_multip(tdm, s, myy, om, 'b');
        re2[2][1]=ex_multip(tdm, s, mzz, om, 'b');
        rhre[0][1]=ex_multip(tdm, mx, mx, om, 'b');
        rhre[1][1]=ex_multip(tdm, my, my, om, 'b');
        rhre[2][1]=ex_multip(tdm, mz, mz, om, 'b');

    }else{

        rh[0][1]=rh[0][0];
        rh[1][1]=rh[1][0];
        rh[2][1]=rh[2][0];
        rh2[0][1]=rh2[0][0];
        rh2[1][1]=rh2[1][0];
        rh2[2][1]=rh2[2][0];
        re[0][1]=re[0][0];
        re[1][1]=re[1][0];
        re[2][1]=re[2][0];
        re2[0][1]=re2[0][0];
        re2[1][1]=re2[1][0];
        re2[2][1]=re2[2][0];
        rhre[0][1]=rhre[0][0];
        rhre[1][1]=rhre[1][0];
        rhre[2][1]=rhre[2][0];
    }//end else

} // end fct

double ex_analyse::ex_d_ex_c (char koord, char spin){
    double erg=0;

    switch (koord){

    case 'x':{
        switch (spin){

        case 'a':{
            erg=sqrt(rh2[0][0]-2*rhre[0][0]+re2[0][0]);
                break;
                }//end case a

        case 'b':{
            erg=sqrt(rh2[0][1]-2*rhre[0][1]+re2[0][1]);
                break;
                }//end case b

        default:break;
        }//end switch

        break;
        }//end case x

    case 'y':{
        switch (spin){

        case 'a':{
            erg=sqrt(rh2[1][0]-2*rhre[1][0]+re2[1][0]);
                break;
                }//end case a

        case 'b':{
            erg=sqrt(rh2[1][1]-2*rhre[1][1]+re2[1][1]);
                break;
                }//end case b

        default:break;
        }//end switch

        break;
        }//end case y

    case 'z':{
        switch (spin){

        case 'a':{
            erg=sqrt(rh2[2][0]-2*rhre[2][0]+re2[2][0]);
                break;
                }//end case a

        case 'b':{
            erg=sqrt(rh2[2][1]-2*rhre[2][1]+re2[2][1]);
                break;
                }//end case b

        default:break;
        }//end switch

        break;
        }//end case z

    default: break;
    }//end switch

    return erg;

}//end fct

double ex_analyse::ex_d_ex_tot (char spin){

    double erg=0;

    switch (spin){

        case 'a':{
            erg=sqrt((ex_d_ex_c('x','a'))*(ex_d_ex_c('x','a'))
                    +(ex_d_ex_c('y','a'))*(ex_d_ex_c('y','a'))
                    +(ex_d_ex_c('z','a'))*(ex_d_ex_c('z','a')));
            break;
        }//end case

        case 'b':{
            erg=sqrt(ex_d_ex_c('x','b')*ex_d_ex_c('x','b')
                    +ex_d_ex_c('y','b')*ex_d_ex_c('y','b')
                    +ex_d_ex_c('z','b')*ex_d_ex_c('z','b'));
            break;
        }//end case

        default:break;
    }//end switch

    return erg;
}//end fct

double ex_analyse::get_rh (char koord, char spin){

    switch (koord){

        case 'x':{
            switch (spin){

            case 'a':{
                return rh[0][0];
                    break;
                    }//end case a

            case 'b':{
                return rh[0][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case x

        case 'y':{
            switch (spin){

            case 'a':{
                return rh[1][0];
                    break;
                    }//end case a

            case 'b':{
                return rh[1][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch
            break;
            }//end case y

        case 'z':{
            switch (spin){

            case 'a':{
                return rh[2][0];
                    break;
                    }//end case a

            case 'b':{
                return rh[2][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case z

        default: return 0;
        }//end switch

}

double ex_analyse::get_re (char koord, char spin){

    switch (koord){

        case 'x':{

            switch (spin){

            case 'a':{
                return re[0][0];
                    break;
                    }//end case a

            case 'b':{
                return re[0][1];
                    break;
                    }//end case b

            default:return 0;

            }//end switch

            break;
            }//end case x

        case 'y':{

            switch (spin){

            case 'a':{
                return re[1][0];
                    break;
                    }//end case a

            case 'b':{
                return re[1][1];
             break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case y

        case 'z':{
            switch (spin){
            case 'a':{
                return re[2][0];
                    break;
                    }//end case a

            case 'b':{
                return re[2][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case z

        default: return 0;
        }//end switch

}

double ex_analyse::get_rh2 (char koord, char spin){

    switch (koord){

        case 'x':{

            switch (spin){

            case 'a':{
                return rh2[0][0];
                    break;
                    }//end case a

            case 'b':{
                return rh2[0][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case x

        case 'y':{

            switch (spin){

            case 'a':{
                return rh2[1][0];
                    break;
                    }//end case a

            case 'b':{
                return rh2[1][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case y

        case 'z':{

            switch (spin){

            case 'a':{
                return rh2[2][0];
                    break;
                    }//end case a

            case 'b':{
                return rh[2][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case z

        default: return 0;
        }//end switch
}

double ex_analyse::get_re2 (char koord, char spin){
    switch (koord){
        case 'x':{
            switch (spin){
            case 'a':{
                return re2[0][0];
                    break;
                    }//end case a
            case 'b':{
                return re2[0][1];
                    break;
                    }//end case b
            default:return 0;
            }//end switch
            break;
            }//end case x
        case 'y':{
            switch (spin){
            case 'a':{
                return re2[1][0];
                    break;
                    }//end case a
            case 'b':{
                return re2[1][1];
                    break;
                    }//end case b
            default:return 0;
            }//end switch
            break;
            }//end case y
        case 'z':{
            switch (spin){
            case 'a':{
                return re2[2][0];
                    break;
                    }//end case a
            case 'b':{
                return re2[2][1];
                    break;
                    }//end case b
            default:return 0;
            }//end switch
            break;
            }//end case z
        default: return 0;
        }//end switch
}

double ex_analyse::get_rhre (char koord, char spin){

    switch (koord){

        case 'x':{
            switch (spin){

            case 'a':{
                return rhre[0][0];
                    break;
                    }//end case a

            case 'b':{
                return rhre[0][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case x

        case 'y':{
            switch (spin){

            case 'a':{
                return rhre[1][0];
                    break;
                    }//end case a

            case 'b':{
                return rhre[1][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case y

        case 'z':{
            switch (spin){

            case 'a':{
                return rhre[2][0];
                    break;
                    }//end case a

            case 'b':{
                return rhre[2][1];
                    break;
                    }//end case b

            default:return 0;
            }//end switch

            break;
            }//end case z

        default: return 0;
        }//end switch
}

double ex_analyse::ex_mean_sep (char spin){

    double vek[3];
    double erg=0;

        switch (spin){

            case 'a':{
                vek[0]=rh[0][0]-re[0][0];
                vek[1]=rh[1][0]-re[1][0];
                vek[2]=rh[2][0]-re[2][0];
                erg=sqrt(vek[0]*vek[0]+vek[1]*vek[1]+vek[2]*vek[2]);
                break;
            }//end case

            case 'b':{
                vek[0]=rh[0][1]-re[0][1];
                vek[1]=rh[1][1]-re[1][1];
                vek[2]=rh[2][1]-re[2][1];
                erg=sqrt(vek[0]*vek[0]+vek[1]*vek[1]+vek[2]*vek[2]);
                break;
            }//end case

            default:break;
        }//end switch

        return erg;
} //end fct

double ex_analyse::ex_sig_h (char spin){

    double erg=0;
    double rh2tot=0;
    double rhtot=0;

        switch (spin){

            case 'a':{
                rh2tot=rh2[0][0]+rh2[1][0]+rh2[2][0];
                rhtot=(rh[0][0])*(rh[0][0])+(rh[1][0])*(rh[1][0])
                        +(rh[2][0])*(rh[2][0]);
                erg=sqrt(rh2tot-rhtot);
                break;
            }//end case

            case 'b':{
                rh2tot=rh2[0][1]+rh2[1][1]+rh2[2][1];
                rhtot=(rh[0][1])*(rh[0][1])+(rh[1][1])*(rh[1][1])
                        +(rh[2][1])*(rh[2][1]);
                erg=sqrt(rh2tot-rhtot);
                break;
            }//end case

            default:break;
            }//end switch

    return erg;
}//end fct

double ex_analyse::ex_sig_e (char spin){
    double erg=0;
    double re2tot=0;
    double retot=0;

        switch (spin){

            case 'a':{
                re2tot=re2[0][0]+re2[1][0]+re2[2][0];
                retot=(re[0][0])*(re[0][0])+(re[1][0])*(re[1][0])
                        +(re[2][0])*(re[2][0]);
                erg=sqrt(re2tot-retot);
                break;
            }//end case

            case 'b':{
                re2tot=re2[0][1]+re2[1][1]+re2[2][1];
                retot=(re[0][1])*(re[0][1])+(re[1][1])*(re[1][1])
                        +(re[2][1])*(re[2][1]);
                erg=sqrt(re2tot-retot);
                break;
            }//end case

            default:break;

            }//end switch

    return erg;

}//end fct

double ex_analyse::ex_cov (char spin){
   double erg=0;

    switch (spin){

            case 'a':{
                erg=rhre[0][0]+rhre[1][0]+rhre[2][0]
                    -(rh[0][0])*(re[0][0])-(rh[1][0])*(re[1][0])
                    -(rh[2][0])*(re[2][0]);
                break;
            }//end case

            case 'b':{
                erg=rhre[0][1]+rhre[1][1]+rhre[2][1]
                    -(rh[0][1])*(re[0][1])-(rh[1][1])*(re[1][1])
                    -(rh[2][1])*(re[2][1]);
                break;
            }//end case

            default:break;
    }//end switch

    return erg;

}//end fct

double ex_analyse::ex_corr(char spin){

    double erg=0;

        switch (spin){

            case 'a':{
                erg=ex_cov('a')/(ex_sig_e('a')*ex_sig_h('a'));
                    break;
                     }//end case

             case 'b':{
                 erg=ex_cov('b')/(ex_sig_e('b')*ex_sig_h('b'));
                    break;
                      }//end case

             default:break;
                   }//end switchfish

        return erg;
}


double ex_analyse::get_sep(char spin) {
    switch (spin) {

    case 'a': {
        return sep[0];
        break;
    } //end case

    case 'b': {
        return sep[1];
        break;
    } //end case

    default:{
        return 0;
        break;
    }//end default
    } //end switchfish
} //end fct

double ex_analyse::get_dex_c(char coord, char spin) {
    switch (coord) {

    case 'x': {
        switch (spin) {

        case 'a': {
            return dex_c[0][0];
            break;
        } //end case a

        case 'b': {
            return dex_c[0][1];
            break;
        } //end case b

        default:
            return 0;
        } //end switch

        break;
    } //end case x

    case 'y': {
        switch (spin) {

        case 'a': {
            return dex_c[1][0];
            break;
        } //end case a

        case 'b': {
            return dex_c[1][1];
            break;
        } //end case b

        default:
            return 0;
        } //end switch

        break;
    } //end case y

    case 'z': {
        switch (spin) {

        case 'a': {
            return dex_c[2][0];
            break;
        } //end case a

        case 'b': {
            return dex_c[2][1];
            break;
        } //end case b

        default:
            return 0;
        } //end switch

        break;
    } //end case z

    default:
        return 0;
    } //end switch
} //end fct

double ex_analyse::get_dex_tot(char spin) {
switch (spin){
case 'a':
    {
        return dex_tot[0];
        break;
    } //end case

    case 'b':
    {
        return dex_tot[1];
        break;
    } //end case

    default:
    {
        return 0;
        break;
    } //end default
} //end switch
} //end fct

double ex_analyse::get_sig_h(char spin) {
    switch(spin){
    case 'a':
    {
        return sig_h[0];
        break;
    } //end case

    case 'b':
    {
        return sig_h[1];
        break;
    } //end case

    default:
    {
        return 0;
        break;
    } //end default
} //end switch
} //end fct

double ex_analyse::get_sig_e(char spin) {
    switch (spin) {
    case 'a': {
        return sig_e[0];
        break;
    } //end case

    case 'b': {
        return sig_e[1];
        break;
    } //end case

    default: {
        return 0;
        break;
    } //end default
    } //end switch
}//end fct

double ex_analyse::get_cov(char spin) {
    switch (spin) {
    case 'a': {
        return cov[0];
        break;
    } //end case

    case 'b': {
        return cov[1];
        break;
    } //end case

    default: {
        return 0;
        break;
    } //end default
    } //end switch
} //end fct

double ex_analyse::get_corr(char spin) {
    switch (spin) {
    case 'a': {
        return corr[0];
        break;
    } //end case

    case 'b': {
        return corr[1];
        break;
    } //end case

    default: {
        return 0;
        break;
    } //end default
    } //end switch
} //end fct

void ex_analyse::perform(const ab_matrix &tdm, const ab_matrix &s,
                const ab_matrix &mxx, const ab_matrix &mx, const ab_matrix &myy,
                const ab_matrix &my, const ab_matrix &mzz,
                const ab_matrix &mz, const ab_matrix &om){
//Forming all needed values for the calculations;
    ex_form(tdm, s, mxx, mx, myy, my, mzz, mz, om);
//Performing all calculations, saving them in the resp. arrays.

    sep[0]=ex_mean_sep('a');
    dex_c[0][0]=ex_d_ex_c ('x', 'a');
    dex_c[1][0]=ex_d_ex_c ('y', 'a');
    dex_c[2][0]=ex_d_ex_c ('z', 'a');
    dex_tot[0]=ex_d_ex_tot('a');
    sig_h[0]=ex_sig_h('a');
    sig_e[0]=ex_sig_e('a');
    cov[0]=ex_cov('a');
    corr[0]=ex_corr('a');

    if (!tdm.is_alpha_eq_beta()){ //If beta is not equal alpha

        sep[1] = ex_mean_sep('b');
        dex_c[0][1] = ex_d_ex_c('x', 'b');
        dex_c[1][1] = ex_d_ex_c('y', 'b');
        dex_c[2][1] = ex_d_ex_c('z', 'b');
        dex_tot[1] = ex_d_ex_tot('b');
        sig_h[1] = ex_sig_h('b');
        sig_e[1] = ex_sig_e('b');
        cov[1] = ex_cov('b');
        corr[1] = ex_corr('b');

    }else{

        sep[1] = sep[0];
        dex_c[0][1] = dex_c[0][0];
        dex_c[1][1] = dex_c[1][0];
        dex_c[2][1] = dex_c[2][0];
        dex_tot[1] = dex_tot[0];
        sig_h[1] = sig_h[0];
        sig_e[1] = sig_e[0];
        cov[1] = cov[0];
        corr[1] = corr[0];

    }//end else
}//end fct



}// end namespace libwfa





