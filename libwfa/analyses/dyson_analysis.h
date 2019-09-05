#ifndef LIBWFA_DYSON_ANALYSIS_H
#define LIBWFA_DYSON_ANALYSIS_H

#include <iostream>
#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/orbital_data.h>
#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {

class orbital_printer_i;

/** \brief Perform Dyson orbital analysis

     \ingroup libwfa
 **/
class dyson_analysis {
private:
    const orbital_data &m_odata; //!< Dyson orbital data

public:
    /** \brief Constructor
        \param odata Dyson orbital data
     **/
    dyson_analysis(const orbital_data &odata) : m_odata(odata) { }

    /** \brief Destructor
     **/
    //~dyson_analysis();

    /** \brief Perform Dyson orbital analysis
        \param out Output stream
        \param ndo # of leading occupation numbers to print (default: 3)
     **/
    void analyse(std::ostream &out, size_t ndo = 3) const;

    /** \brief Export Dyson orbitals
        \param pr Orbital printer
        \param thresh Threshold to determine important orbitals

        Only Dyson orbitals s with occupation number larger than the threshold are exported
     **/
    void export_orbitals(orbital_printer_i &pr, double thresh) const;

private:
    static void build_selector(const arma::vec &e, double thresh,
            orbital_selector &sel);
};


} // namespace libwfa

#endif // LIBWFA_DYSON_ANALYSIS_H
