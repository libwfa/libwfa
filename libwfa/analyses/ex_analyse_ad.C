#include "ex_analyse_ad.h"

namespace libwfa {

using namespace arma;


void ex_analyse_ad::ex_form_ad(const ab_matrix &att, const ab_matrix &det,
        const Mat<double> &s, const contract_ad_i &name){

    prom[0]=trace(att.alpha()*s);
    prom[1]=trace(att.beta()*s);

    rh[0][0]=name.perform(det, prom[0],"x",'a');
    rh[1][0]=name.perform(det, prom[0],"y",'a');
    rh[2][0]=name.perform(det, prom[0],"z",'a');
    rh2[0][0]=name.perform(det, prom[0],"xx",'a');
    rh2[1][0]=name.perform(det, prom[0],"yy",'a');
    rh2[2][0]=name.perform(det, prom[0],"zz",'a');
    re[0][0]=name.perform(att, prom[0],"x",'a');
    re[1][0]=name.perform(att, prom[0],"y",'a');
    re[2][0]=name.perform(att, prom[0],"z",'a');
    re2[0][0]=name.perform(att, prom[0],"xx",'a');
    re2[1][0]=name.perform(att, prom[0],"yy",'a');
    re2[2][0]=name.perform(att, prom[0],"zz",'a');

    if (!att.is_alpha_eq_beta()){

        rh[0][1]=name.perform(det, prom[1],"x",'b');
        rh[1][1]=name.perform(det, prom[1],"y",'b');
        rh[2][1]=name.perform(det, prom[1],"z",'b');
        rh2[0][1]=name.perform(det, prom[1],"xx",'b');
        rh2[1][1]=name.perform(det, prom[1],"yy",'b');
        rh2[2][1]=name.perform(det, prom[1],"zz",'b');
        re[0][1]=name.perform(att, prom[1],"x",'b');
        re[1][1]=name.perform(att, prom[1],"y",'b');
        re[2][1]=name.perform(att, prom[1],"z",'b');
        re2[0][1]=name.perform(att, prom[1],"xx",'b');
        re2[1][1]=name.perform(att, prom[1],"yy",'b');
        re2[2][1]=name.perform(att, prom[1],"zz",'b');

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
        const Mat<double> &s, const contract_ad_i &name){
//Forming all needed values for the calculations;
    ex_form_ad(att,det,s,name);
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
