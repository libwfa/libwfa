#ifndef EX_ANA_PRINTER_AD_H_
#define EX_ANA_PRINTER_AD_H_

#include <cstddef>
#include <string>
#include <libwfa/analyses/ex_analyse_ad.h>


namespace libwfa{

using namespace arma;

class ex_ana_printer_ad{
private:
    static const double k_au2ang;

    /**\brief Prints the results of the exciton analysis.
        \param[in] spin Print specific spin
        \param[in] analyse_ad Name of the class ex_analyse_ad
        \param[in] out Ostream name
     **/
    void do_print (char spin, ex_analyse_ad &analyse_ad, std::ostream &out);

public:
    /** \brief Calls the print function.
        \param[in] analyse_ad Name of the class ex_analyse_ad
        \param[in] out Ostream name
     **/
    void perform (ex_analyse_ad &analyse_ad, std::ostream &out);

};//end class


}//end namespace libwfa




#endif /* EX_ANA_PRINTER_AD_H_ */
