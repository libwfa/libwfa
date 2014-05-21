#ifndef LIBWFA_CTNUM_EXPORT_H
#define LIBWFA_CTNUM_EXPORT_H

#include <cstddef>
#include <string>
#include "ctnum_printer_i.h"

namespace libwfa {

/** \brief Printer for CT number data (export to file)

    This printer exports the CT number data obtained from one transition density
    to one or two files in the current directory.

    \ingroup libwfa
 **/
class ctnum_export : public ctnum_printer_i {
protected:    
    std::ostream &m_out; //!< Output stream
    size_t m_ncols; //!< Max number of columns
    size_t m_colwidth; //!< Column width
    size_t m_prec; //!< Precision of printed numbers

    std::string m_prefix; //!< File prefix for the data
    double m_energy; //!< Excitation energy of the state
    double m_osc_strength; //!< Oscillator strength of the state

public:
    /** \brief Constructor
        \param ncols Max number of columns
        \param colwidth Max column width
        \param prec Precision the data columns
     **/
    ctnum_export(std::ostream &out, const std::string &prefix = "ctnum",
        size_t ncols = 3, size_t colwidth = 15, size_t prec = 6):
        m_out(out), m_ncols(ncols), m_colwidth(colwidth), m_prec(prec),
        m_prefix(prefix), m_energy(0.0), m_osc_strength(0.0) { }

    /** \brief Set prefix
     **/
    void set_prefix(const std::string &prefix) {
        m_prefix = prefix;
    }

    /** \brief Set state information
     **/
    void set_state_info(double energy, double osc) {
        m_energy = energy; m_osc_strength = osc;
    }

    /** \brief Export the CT number data

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
    virtual void perform(const ab_matrix &om, const double (&om_tot)[2]);

private:
    void do_export(const std::string &fname, const arma::Mat<double> &ct);
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_EXPORT_H
