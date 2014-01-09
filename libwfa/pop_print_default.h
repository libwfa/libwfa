#ifndef LIBWFA_POP_PRINT_DEFAULT_H
#define LIBWFA_POP_PRINT_DEFAULT_H

#include <cstddef>
#include <iostream>
#include "pop_print_i.h"

namespace libwfa {


/** \brief Printer of N sets of population data

 **/
class pop_print_default : public pop_print_i {
private:
    const std::vector<std::string> &m_labels; //!< Labels for each line
    std::ostream &m_out; //!< Output stream
    size_t m_colwidth; //!< Column width
    size_t m_prec; //!< Precision of printed numbers
    
public:
    pop_print_default(const std::vector<std::string> &l,
        std::ostream &out, size_t colwidth = 20, size_t prec = 6) :
        m_labels(l), m_out(out), m_colwidth(colwidth), m_prec(prec) { }

    virtual void perform(const pop_data &p);
};


} // namespace libwfa

#endif // LIBWFA_POP_PRINT_DEFAULT_H
