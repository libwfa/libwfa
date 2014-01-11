#ifndef LIBWFA_CTNUM_EXPORT_H
#define LIBWFA_CTNUM_EXPORT_H

#include <cstddef>
#include <string>
#include "ctnum_print_i.h"

namespace libwfa {

/** \brief Printer for several sets of CT number data (export to file)

    This printer exports each data set to a separate file in the current
    directory.

    \ingroup libwfa
 **/
class ctnum_export : public ctnum_print_i {
protected:    
    std::string m_prefix; //!< Prefix for each file
    size_t m_ncols; //!< Max number of columns
    size_t m_colwidth; //!< Column width
    size_t m_prec; //!< Precision of printed numbers

public:
    /** \brief Constructor
        \param prefix Prefix attached to each file name
        \param ncols Max number of columns
        \param colwidth Max column width
        \param prec Precision the data columns
     **/
    ctnum_export(const std::string &prefix = "",
        size_t ncols = 3, size_t colwidth = 15, size_t prec = 6):
        m_prefix(prefix), m_ncols(ncols), m_colwidth(colwidth), m_prec(prec) { }
    
    /** \brief Export the CT number data

        The CT number data is exported as one text file per data set into the
        current directory. The files are named as <prefix><name>.om, thus the
        user should ensure that no two data sets have the same name.
        The format of the text files is as follows:
        - the first line contains name, excitation energy and oscillator
          strength each separated by spaces
        - the second line contains the dimensions of the data (always 2) and
          the number of rows and columns of the data
        - in the remaining lines the CT number data is printed as linear
          array with fixed precision numbers using the provided width and
          precision. A line break is inserted after ncols numbers have been
          written.

        The exported data can be analysed further using the python script
        \c run_dens_ana.py available at
        http://www.iwr.uni-heidelberg.de/groups/compchem/personal/felix_plasser/download.html
     **/
    virtual void perform(const ctnum_data &ct);

};

} // namespace libwfa

#endif // LIBWFA_CTNUM_EXPORT_H
