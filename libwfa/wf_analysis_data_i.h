#ifndef LIBWFA_WF_ANALYSIS_DATA_I_H
#define LIBWFA_WF_ANALYSIS_DATA_I_H

#include <libwfa/analyses/ctnum_analysis_i.h>
#include <libwfa/analyses/pop_analysis_i.h>
#include <libwfa/core/mom_builder_i.h>
#include <libwfa/export/ctnum_printer_i.h>
#include <libwfa/export/density_printer_i.h>
#include <libwfa/export/orbital_printer_i.h>
#include <libwfa/soc.h>

namespace libwfa {


/** \brief Interface to provide the data required to perform the analysis

    If new analyses are added to libwfa which require additional data, please
    add the respective functions to the interface.

    \ingroup libwfa
 **/
class wf_analysis_data_i {
public:
    /** Types of wave function / density analyses to be performed
     **/
    // PP: Dear developers,
    // If you want to add a new analysis_type, 
    // please check that the number of analysis types corresponds 
    // to all the available types. Otherwise, this will happen:
    // terminate called after throwing an instance of 'std::out_of_range'
    //   what():  bitset::set: __position (which is 10) >= _Nb (which is 10)
    // Exception will be thrown in 
    // qchem/qchem_wf_analysis_data.C
    // See #1797 ticket on trac as an example of how it looks like
    // Thank you.
    enum analysis_type {
        NO = 0,     //!< Natural orbital analysis
        DENS_MOM,   //!< Multipole moments of the density
        NDO,        //!< Natural difference orbital analysis
        FORM_AD,    //!< Form attachment/detachment densities
        EXCITON_AD, //!< Exciton analysis on a/d densities
        NTO,        //!< Natural transition orbital analysis
        TDEN_MOM,   //!< Multipole moments of the transition density
        FORM_EH,    //!< Form electron/hole densities
        EXCITON,    //!< Exciton analysis on transition density
        SA_NTO,     //!< State-averaged NTO analysis
        NTO_ENE,    //!< Compute "energies" of hole and particle NTOs
        DYSON,      //!< Dyson orbital analysis
	NTO_SOC     //!< Compute SOC integrals of hole and particle NTOs 
                    //!<  of the universal triplet optdm
    };
    enum { WFA_TYPES = 13 }; //!< Number of analysis types

    /** \brief Orbital parameters
     **/
    struct orbital_params {
        size_t norb; //!< Number of leading orbitals
        double thresh; //!< Threshold of important orbitals

        /** \brief Constructor
         **/
        orbital_params(size_t norb_ = 0, double thresh_ = 0.0) :
            norb(norb_), thresh(thresh_) { }
    };

public:
    /** \brief Virtual destructor
     **/
    virtual ~wf_analysis_data_i() { }

    /** \brief Return if a certain analysis should be performed
     **/
    virtual bool is_active(enum analysis_type t) = 0;

    /** \brief Return parameters for orbital analyses
        \return Pair comprising the number of leading orbitals to print and
            an threshold for important orbitals (see e.g. \ref no_analysis
            for details)
     **/
    virtual orbital_params get_orbital_params(enum orbital_type::ot t) = 0;

    /** \brief Retrieve the AO overlap matrix
     **/
    virtual const arma::mat &overlap() = 0;

    /** \brief Retrieve the MO coefficient matrices
     **/
    virtual const ab_matrix &coefficients() = 0;

    /** \brief Fock matrix in AO
     **/
    virtual const ab_matrix &fock() = 0;

    /** \brief 1e SOC operators in AO
     **/
    virtual const h_so &so1e() = 0;

    /** \brief mean-field SOC operators in AO
     **/
    virtual const h_so &somf() = 0;

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
    virtual std::unique_ptr<ctnum_printer_i> ctnum_printer(size_t i,
            const std::string &name, const std::string &desc) = 0;

    /** \brief Name of i-th CT number analysis
     **/
    virtual const std::vector<std::string> &prop_list() = 0;

    //@}

    /** \brief Builder of exciton moments
     **/
    virtual const mom_builder_i &mom_builder() = 0;

    /** \brief Return the atomic coordinates
     **/
    virtual const arma::mat &coordinates() = 0;

    /** \brief Return the atomic numbers
     **/
    virtual const arma::vec &atomic_charges() = 0;
};


} // namespace libwfa

#endif // LIBWFA_WF_ANALYSIS_DATA_I_H
