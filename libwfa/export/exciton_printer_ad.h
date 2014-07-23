#ifndef LIBWFA_EXCITON_PRINTER_AD_H
#define LIBWFA_EXCITON_PRINTER_AD_H

#include <libwfa/analyses/exciton_moments.h>
#include <ostream>
#include "exciton_printer_i.h"

namespace libwfa{

/** \brief Printer of exciton moment data (A/D)

    \ingroup libwfa
 **/
class exciton_printer_ad : public exciton_printer_i {
public:
    /** \copydoc exciton_printer_i::perform
     **/
    void perform(ab_exciton_moments &mom, std::ostream &out) const;

private:
    /** \brief Prints the exciton moment data for one spin
        \param mom Moment data
        \param out Output stream
     **/
    void print(exciton_moments &mom, std::ostream &out) const;
};


}//end namespace libwfa




#endif // LIBWFA_EXCITON_PRINTER_AD_H
