#ifndef LIBWFA_CTNUM_EXPORT_H
#define LIBWFA_CTNUM_EXPORT_H

#include <cstddef>
#include <string>
#include <libwfa/core/ab_matrix.h>
#include "ctnum_printer_i.h"

namespace libwfa {

/** \brief Printer for CT number data (export to file)

    This printer exports the CT number data obtained from one transition density
    to one or two files in the current directory.

    \ingroup libwfa
 **/
class ctnum_export : public ctnum_printer_i {
private:
    size_t m_ncols; //!< Max number of columns
    size_t m_colwidth; //!< Column width
    size_t m_prec; //!< Precision of printed numbers

    std::string m_sid; //!< State ID (used as file prefix)
    std::string m_sdesc; //!< State description
    double m_energy; //!< Excitation energy of the state
    double m_osc_strength; //!< Oscillator strength of the state

public:
    /** \brief Constructor
        \param ncols Max number of columns
        \param colwidth Max column width
        \param prec Precision the data columns
     **/
    ctnum_export(size_t ncols = 3, size_t colwidth = 15, size_t prec = 6):
        m_ncols(ncols), m_colwidth(colwidth), m_prec(prec),
        m_energy(0.0), m_osc_strength(0.0) { }

    /** \brief Set state information
        \param sid State ID
        \param sdesc State descripion
        \param energy State energy
        \param osc Oscillator strength
     **/
    void set_state_info(const std::string &sid,
        const std::string &sdesc, double energy, double osc);

    /** \brief Export the CT number data
        \param om CT number data (omega matrix)
        \param om_tot Omega total
        \param out Output stream

        The CT number data which is passed as argument ct is exported into one
        or two files. If ct contains one matrix (alpha == beta), one file is
        created with the name \<prefix\>.om. Otherwise the data is exported
        as two files named \<prefix\>(_a|_b).om. To avoid overwriting files
        the user should ensure that the prefix is changed after every call of
        perform. The format of the text files is as follows:
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
    virtual void perform(const ab_matrix &om,
        const double (&om_tot)[2], std::ostream &out) const;

private:
    void do_export(const std::string &fname, const arma::Mat<double> &ct) const;
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_EXPORT_H
