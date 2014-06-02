#ifndef EX_ANA_PRINTER_H_
#define EX_ANA_PRINTER_H_

#include <cstddef>
#include <string>
#include <libwfa/analyses/ex_analyse.h>


namespace libwfa{

using namespace arma;

class ex_ana_printer{
private:
    /**\brief Prints the results of the exciton analysis.
         \param[in] aeqb Alph==beta?
             \param[in] analyse Name of the class ex_analyse
             \param[in] out Ostream name
     **/

    void do_print (bool aeqb, ex_analyse &analyse, std::ostream &out);
public:
    /** \brief Calls the print function.
             \param[in] aeqb Alph==beta?
             \param[in] analyse Name of the class ex_analyse
             \param[in] out Ostream name
     **/
    void perform (bool aeqb, ex_analyse &analyse, std::ostream &out);

};//end class


}//end namespace libwfa



#endif /* EX_ANA_PRINTER_H_ */
