#ifndef LIBWFA_POP_PRINT_I_H
#define LIBWFA_POP_PRINT_I_H

#include <list>
#include <string>
#include <vector>
#include <utility>

namespace libwfa {


/** \brief Interface for printer of population data

    \ingroup libwfa
 **/
class pop_print_i {
public:
    typedef std::list< std::pair<std::string, std::vector<double> > > pop_data;

public:
    virtual ~pop_print_i() { }

    virtual void perform(pop_data &p) = 0;
};


} // namespace adcman

#endif // LIBWFA_POP_PRINT_I_H
