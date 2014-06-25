#ifndef LIBWFA_POP_PRINTER_DEFAULT_H
#define LIBWFA_POP_PRINTER_DEFAULT_H

#include <cstddef>
#include <iostream>
#include "pop_printer_i.h"

namespace libwfa {


/** \brief Printer of N sets of population data to output stream

    The population data is printed as one table

    \ingroup libwfa
 **/
class pop_printer_default : public pop_printer_i {
private:
    const std::vector<std::string> &m_labels; //!< Labels for each line
    size_t m_colwidth; //!< Column width
    size_t m_prec; //!< Precision of printed numbers
    size_t m_offset; //!< Offset of table
    
public:
    /** \brief Constructor
        \param l Labels for each line of the table (= length of each data set
        \param out Output stream
        \param colwidth Max column width
        \param prec Precision the data columns
     **/
    pop_printer_default(const std::vector<std::string> &l,
        size_t colwidth = 16, size_t prec = 6, size_t offset = 2) :
        m_labels(l), m_colwidth(colwidth), m_prec(prec), m_offset(offset) { }

    /** \brief Print the population data

        Population data is printed as one table and each number with the
        provided precision. The method attempts to fit the table into
        lines of 80 character by decreasing the provided column width,
        if necessary and possible (due to precision). However, it does not
        introduce line breaks, if the table exceeds the total width of 80
        characters.
     **/
    virtual void perform(const pop_data &p, std::ostream &out) const;
};


} // namespace libwfa

#endif // LIBWFA_POP_PRINTER_DEFAULT_H
