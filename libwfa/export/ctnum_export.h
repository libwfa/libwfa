#ifndef LIBWFA_CTNUM_EXPORT_H
#define LIBWFA_CTNUM_EXPORT_H

#include <cstddef>
#include <string>
#include <libwfa/core/ab_matrix.h>
#include "ctnum_printer_i.h"

namespace libwfa {

/** \brief Printer for CT number data (export to file)

    This printer exports the CT number data obtained from one transition
    density to one or two files in the current directory.

    \ingroup libwfa
 **/
class ctnum_export : public ctnum_printer_i {
private:
    std::string m_prefix; //!< File prefix
    std::string m_desc; //!< State description

    size_t m_ncols; //!< Max number of columns
    size_t m_colwidth; //!< Column width
    size_t m_prec; //!< Precision of printed numbers

public:
    /** \brief Constructor
        \param prefix File prefix
        \param desc Description
        \param ncols Max number of columns
        \param colwidth Max column width
        \param prec Precision the data columns
     **/
    ctnum_export(const std::string &prefix, const std::string &desc,
        size_t ncols = 3, size_t colwidth = 15, size_t prec = 6) :
        m_prefix(prefix), m_desc(desc),
        m_ncols(ncols), m_colwidth(colwidth), m_prec(prec) { }

    /** \brief Export the CT number data
        \param om CT number data (omega matrix)

        The CT number data which is passed as argument ct is exported into one
        or two files. If ct contains one matrix (alpha == beta), one file is
        created with the name \<prefix\>.om. Otherwise the data is exported
        as two files named \<prefix\>(_a|_b).om. To avoid overwriting files
        the user should ensure that the prefix is changed after every call of
        perform. The format of the text files is as follows:
        - the first line contains the description which should be composed of
          name, energy, and oscillator strength of the state
        - the second line contains the dimensions of the data (always 2) and
          the number of rows and columns of the data
        - in the remaining lines the CT number data is printed as linear
          array with fixed precision numbers using the provided width and
          precision. A line break is inserted after ncols numbers have been
          written.

        The exported data can be analysed further using TheoDORE
        available at
        http://theodore-qc.sourceforge.net/
     **/
    virtual void perform(const ab_matrix &om) const;

private:
    void do_export(const std::string &fname, const arma::mat &ct) const;
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_EXPORT_H
