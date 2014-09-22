#ifndef LIBWFA_WF_ANALYSIS_DATA_I_H
#define LIBWFA_WF_ANALYSIS_DATA_I_H

#include <libwfa/analyses/ctnum_analysis_i.h>
#include <libwfa/analyses/pop_analysis_i.h>
#include <libwfa/core/mom_builder_i.h>
#include <libwfa/export/ctnum_printer_i.h>
#include <libwfa/export/density_printer_i.h>
#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {


/** \brief Interface to provide the data required to perform the analysis

    If new analyses are added to libwfa which require additional data, please
    add the respective functions to the interface.

    \ingroup libwfa
 **/
class wf_analysis_data_i {
public:
    /** \brief Virtual destructor
     **/
    virtual ~wf_analysis_data_i() { }

    /** \brief Retrieve the AO overlap matrix
     **/
    virtual const arma::mat &overlap() = 0;

    /** \brief Retrieve the MO coefficient matrices
     **/
    virtual const ab_matrix &coefficients() = 0;

    /** \brief Construct a printer of density matrices
        \param name Name of state to which the density matrices belong
            (should be usable as file name)
        \param desc Description of state
        \return Pointer to new density printer
     **/
    virtual density_printer_i *density_printer(const std::string &name,
            const std::string &desc) = 0;

    /** \brief Construct a printer of orbitals
        \param name Name of state to which the orbitals belong
            (should be usable as file name)
        \param desc Description of state
        \return Pointer to new orbital printer
     **/
    virtual orbital_printer_i *orbital_printer(const std::string &name,
            const std::string &desc) = 0;

    //! \name Population analysis related functions
    //@{

    /** \brief Number of population analyses available
     **/
    virtual size_t n_pop_analyses() = 0;

    /** \brief Name of i-th population analysis
     **/
    virtual const std::string &pop_name(size_t i) = 0;

    /** \brief Row labels for i-th population analysis
     **/
    virtual const std::vector<std::string> &pop_labels(size_t i) = 0;

    /** \brief i-th population analysis
     **/
    virtual const pop_analysis_i &pop_analysis(size_t i) = 0;

    /** \brief Base populations for i-th population analysis
     **/
    virtual const arma::vec &ref_population(size_t i) = 0;

    //@}

    //! \name CT number analysis related functions
    //@{

    /** \brief Number of CT number analyses available
     **/
    virtual size_t n_ctnum_analyses() = 0;

    /** \brief Name of i-th CT number analysis
     **/
    virtual const std::string &ctnum_name(size_t i) = 0;

    /** \brief i-th CT number analysis
     **/
    virtual const ctnum_analysis_i &ctnum_analysis(size_t i) = 0;

    /** \brief Printer of i-th CT number data
     **/
    virtual std::auto_ptr<ctnum_printer_i> ctnum_printer(size_t i,
            const std::string &name, const std::string &desc) = 0;

    //@}

    /** \brief Builder of exciton moments
     **/
    virtual const mom_builder_i &mom_builder() = 0;
};


} // namespace libwfa

#endif // LIBWFA_WF_ANALYSIS_DATA_I_H
