#ifndef EX_ANA_PRINTER_H_
#define EX_ANA_PRINTER_H_

#include <cstddef>
#include <string>
#include <libwfa/analyses/ex_analyse.h>


namespace libwfa{

using namespace arma;

class ex_ana_printer{
private:
	static const double k_au2ang;

public:
    /** \brief Calls the print function.
        \param[in] analyse Analysis class
        \param[in] out Ostream name
     **/
    void perform(ex_analyse &analyse, std::ostream &out);

private:
    /** \brief Calls the print function.
    	\param[in] spin Print specific spin
        \param[in] analyse Analysis class
        \param[in] out Ostream name
     **/
    void do_print(char spin, ex_analyse &analyse, std::ostream &out);
};//end class


}//end namespace libwfa



#endif /* EX_ANA_PRINTER_H_ */
