#ifndef LIBWFA_WF_ANALYSIS_DATA_I_H
#define LIBWFA_WF_ANALYSIS_DATA_I_H

#include <libwfa/analyses/ctnum_analysis_i.h>
#include <libwfa/analyses/pop_analysis_i.h>
#include <libwfa/core/mom_builder_i.h>
#include <libwfa/export/ctnum_printer_i.h>
#include <libwfa/export/density_printer_i.h>
#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {


/** \brief Parameters of state density matrix analysis
 **/
struct opdm_params {
    bool no_analysis; //!< Perform NO analysis
    size_t nno; //!< Number of leading NOs to print
};


/** \brief Parameters of difference density matrix analysis
 **/
struct opddm_params {
    /** \brief Parameter for NDO analysis
     **/
    typedef enum {
        NONE = 0, //!< Do nothing
        NDO = 1, //!< Perform NDO analysis only
        FORM_AD = 2, //!< Form attachment / detachment densities
        EXCITON = 3 //!< Perform exciton analysis on top of A/D densities
    } ndo_param;

    ndo_param ndo_analysis; //!< NDO analysis parameter
    size_t nndo; //!< Number of leading NDO pairs to print
};


/** \brief Parameters of transition density matrix analysis
 **/
struct optdm_params {
    /** \brief Parameter for electron / hole density analysis
     **/
    typedef enum {
        NONE = 0, //!< Do nothing
        FORM_EH = 1, //!< Form electron / hole densities
        AVERAGE = 2 //!< Add to average electron / hole densities
    } eh_param;

    bool nto_analysis; //!< Perform NTO analysis
    size_t nnto; //!< Number of leading NTOs to print
    double thresh; //!< Zero threshold for NTO occupation numbers
    eh_param eh; //!< Electron / hole densities parameter
    bool exciton_analysis; //!< Perform exciton analysis
};


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
    virtual std::auto_ptr<density_printer_i> density_printer(
            const std::string &name, const std::string &desc) = 0;

    /** \brief Construct a printer of orbitals
        \param name Name of state to which the orbitals belong
            (should be usable as file name)
        \param desc Description of state
        \return Pointer to new orbital printer
     **/
    virtual std::auto_ptr<orbital_printer_i> orbital_printer(
            const std::string &name, const std::string &desc) = 0;

    /** \brief Parameter structure for state density matrix analyses
     **/
    virtual opdm_params  get_opdm_params() = 0;

    /** \brief Parameter structure for difference density matrix analyses
     **/
    virtual opddm_params get_opddm_params() = 0;

    /** \brief Parameter structure for transition density matrix analyses
     **/
    virtual optdm_params get_optdm_params() = 0;

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
