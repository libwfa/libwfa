#ifndef LIBWFA_EXCITON_PRINTER_I_H
#define LIBWFA_EXCITON_PRINTER_I_H

#include <libwfa/analyses/exciton_moments.h>
#include <ostream>

namespace libwfa{

class ab_exciton_moments;

/** \brief Interface for printer of exciton moments / information

    \ingroup libwfa
 **/
class exciton_printer_i{
public:
    /** \brief Virtual destructor
     **/
    virtual ~exciton_printer_i() { }

    /** \brief Prints the exciton moments
        \param mom Exciton moment data
        \param out Output stream
     **/
    virtual void perform(ab_exciton_moments &mom, std::ostream &out) const = 0;
};

}//end namespace libwfa

#endif // LIBWFA_EXCITON_PRINTER_I_H
