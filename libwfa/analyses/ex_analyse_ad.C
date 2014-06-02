#include "ex_analyse_ad.h"

namespace libwfa {

using namespace arma;

double ex_analyse_ad::ex_multip_ad (const ab_matrix &ad,
        const ab_matrix &op, const ab_matrix &om, char spin){

    switch (spin){
    case 'a':{
            const Mat<double> &ad_a=ad.alpha();
            const Mat<double> &op_a=op.alpha();
            double omval_a=accu(om.alpha());
            double expval_a=0;
            Mat<double> exp_a=ad_a*op_a;
            expval_a=trace(exp_a)/omval_a;
            return expval_a;
            break;
     }

    case 'b':{
        const Mat<double> &ad_b=ad.beta();
        const Mat<double> &op_b=op.beta();
        double omval_b=accu(om.beta());
        double expval_b=0;
        Mat<double> exp_b=ad_b*op_b;
        expval_b=trace(exp_b)/omval_b;
        return expval_b;
        break;
    }
        default:
            return 0;

    } //end switch

}//end fct

double ex_analyse_ad::ex_form_ad(const ab_matrix &att, const ab_matrix &det,
        const ab_matrix &mx, const ab_matrix &mxx, const ab_matrix &my, const ab_matrix &myy,
        const ab_matrix &mz, const ab_matrix &mzz, const ab_matrix &om){

    rh[0][0]=ex_multip_ad(det, mx, om, 'a');
    rh[1][0]=ex_multip_ad(det, my, om, 'a');
    rh[2][0]=ex_multip_ad(det, mz, om, 'a');
    rh2[0][0]=ex_multip_ad(det, mxx, om, 'a');
    rh2[1][0]=ex_multip_ad(det, myy, om, 'a');
    rh2[2][0]=ex_multip_ad(det, mzz, om, 'a');
    re[0][0]=ex_multip_ad(att, mx, om, 'a');
    re[1][0]=ex_multip_ad(att, my, om, 'a');
    re[2][0]=ex_multip_ad(att, mz, om, 'a');
    re2[0][0]=ex_multip_ad(att, mxx, om, 'a');
    re2[1][0]=ex_multip_ad(att, myy, om, 'a');
    re2[2][0]=ex_multip_ad(att, mzz, om, 'a');

    if (!att.is_alpha_eq_beta()){

        rh[0][1]=ex_multip_ad(det, mx, om, 'b');
        rh[1][1]=ex_multip_ad(det, my, om, 'b');
        rh[2][1]=ex_multip_ad(det, mz, om, 'b');
        rh2[0][1]=ex_multip_ad(det, mxx, om, 'b');
        rh2[1][1]=ex_multip_ad(det, myy, om, 'b');
        rh2[2][1]=ex_multip_ad(det, mzz, om, 'b');
        re[0][1]=ex_multip_ad(att, mx, om, 'b');
        re[1][1]=ex_multip_ad(att, my, om, 'b');
        re[2][1]=ex_multip_ad(att, mz, om, 'b');
        re2[0][1]=ex_multip_ad(att, mxx, om, 'b');
        re2[1][1]=ex_multip_ad(att, myy, om, 'b');
        re2[2][1]=ex_multip_ad(att, mzz, om, 'b');

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
    }//end else


}//end fct

double ex_analyse_ad::ex_mean_sep_ad (char spin){

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

double ex_analyse_ad::ex_sig_h_ad (char spin){

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

double ex_analyse_ad::ex_sig_e_ad (char spin){
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

void ex_analyse_ad::perform(const ab_matrix &att, const ab_matrix &det,
        const ab_matrix &mx,const ab_matrix &mxx, const ab_matrix &my,
        const ab_matrix &myy, const ab_matrix &mz, const ab_matrix &mzz,
        const ab_matrix &om){
//Forming all needed values for the calculations;
    ex_form_ad(att, det, mx, mxx, my, myy, mz, mzz, om);
//Performing all calculations, saving them in the resp. arrays.

    sep[0]=ex_mean_sep_ad('a');
    sig_h[0]=ex_sig_h_ad('a');
    sig_e[0]=ex_sig_e_ad('a');


    if (!att.is_alpha_eq_beta()){ //If beta is not equal alpha

        sep[1] = ex_mean_sep_ad('b');
        sig_h[1] = ex_sig_h_ad('b');
        sig_e[1] = ex_sig_e_ad('b');

    }else{

        sep[1] = sep[0];
        sig_h[1] = sig_h[0];
        sig_e[1] = sig_e[0];

    }//end else
}//end fct

double ex_analyse_ad::get_sig_h(char spin) {
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

double ex_analyse_ad::get_sig_e(char spin) {
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

double ex_analyse_ad::get_sep(char spin) {
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

double ex_analyse_ad::get_rh (char koord, char spin){

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

}//end fct

double ex_analyse_ad::get_re (char koord, char spin){

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

}//end fct

double ex_analyse_ad::get_rh2 (char koord, char spin){

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
}//end fct

double ex_analyse_ad::get_re2 (char koord, char spin){
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
}//end fct


} //end namespace
