#ifndef LIBWFA_NO_ANALYSIS_H
#define LIBWFA_NO_ANALYSIS_H

#include <iostream>
#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/orbital_data.h>
#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {

class orbital_printer_i;

/** \brief Perform NO analysis of a state density matrix

     \ingroup libwfa
 **/
class no_analysis {
private:
    orbital_data *m_no[3]; //!< NO orbitals

public:
    /** \brief Constructor
        \param s Overlap matrix
        \param c Arbitrary coefficient matrix to transform in orthogonal basis
        \param sdm State density matrix
     **/
    no_analysis(const arma::mat &s, const ab_matrix &c, const ab_matrix &sdm);

    /** \brief Destructor
     **/
    ~no_analysis();

    /** \brief Perform NO analysis
        \param out Output stream
        \param nno # of leading occupation numbers to print (default: 3)
     **/
    void analyse(std::ostream &out, size_t nno = 3) const;

    /** \brief Export NOs
        \param pr Orbital printer
     **/
    void export_orbitals(orbital_printer_i &pr) const;

private:
    static void analysis_p1(std::ostream &out,
            const arma::vec &e, size_t nno = 3);

    static void analysis_p2(std::ostream &out, const arma::vec &e);

    static void build_selector(const arma::vec &e, orbital_selector &sel);
};


} // namespace libwfa

#endif // LIBWFA_NO_ANALYSIS_H
