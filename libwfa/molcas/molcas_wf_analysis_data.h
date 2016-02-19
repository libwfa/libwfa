#ifndef LIBWFA_MOLCAS_WF_ANALYSIS_DATA_H
#define LIBWFA_MOLCAS_WF_ANALYSIS_DATA_H

#include "H5Cpp.h"
#include <libwfa/wf_analysis_data_i.h>
#include "molcas_mom_builder.h"

namespace libwfa {
    
class molcas_wf_analysis_data : public libwfa::wf_analysis_data_i {
private:
    /** \brief Internal structure to hold population analyses
     **/
    struct pa_data {
        std::string name;
        const std::vector<std::string> &labels;
        const arma::vec &p0;
        std::auto_ptr<libwfa::pop_analysis_i> analysis;

        pa_data(const std::string &n, const std::vector<std::string> &l,
            libwfa::pop_analysis_i *a, const arma::vec &p) :
            name(n), labels(l), analysis(a), p0(p) { }
     };

    /** \brief Internal structure to hold CT number analyses
     **/
    struct cta_data {
        std::string name;
        std::string suffix;
        std::auto_ptr<libwfa::ctnum_analysis_i> analysis;

        cta_data(const std::string &n, const std::string &s,
            libwfa::ctnum_analysis_i *a) :
            name(n), suffix(s), analysis(a) { }
    };

    /** \brief Internal structure to hold basic molecular data
     **/
    struct base_data {
        std::vector<std::string> atoms; //!< Atom names
        arma::uvec atomic_numbers; //!< Atomic numbers
        arma::vec atomic_charges; //!< Atomic charges
        arma::mat coordinates; //!< Atom coordinates
        arma::uvec bf2atoms; //!< Map of basis functions to atoms
        ab_matrix c_fb; //!< MO coefficients
        arma::mat s; //!< AO-Overlap matrix
        molcas_mom_builder mom; //!< Moments builder
        std::string mo_types_a; //!< Alpha types: F(rozen), I(nactive), (RAS)1,2,3, S(econdary)
        std::string mo_types_b; //!< Beta types: F(rozen), I(nactive), (RAS)1,2,3, S(econdary)

        /** \brief Constructor
            \param nao Number of atomic basis functions
            \param nb2 Number of shell pairs
            \param maxmm Maximum multipole moment available
            \param unr Unrestricted calculation
         **/
        base_data(size_t nao, size_t nb2, size_t maxmm, bool unr) :
            c_fb(unr), mom(nao, nb2, maxmm) { }
    };

    /** \brief Export types
     **/
    enum export_type {
        EXPORT_NONE = 0,//!< EXPORT_NONE
        EXPORT_H5,  //!< EXPORT_H5
    };

public:
    static const char k_clazz[]; //!< Class name

private:
    H5::H5File m_file; //!< HDF5 file
    std::vector<pa_data *> m_pa; //!< Population analyses
//    std::vector<cta_data *> m_cta; //!< CT number analyses

    std::bitset<WFA_TYPES> m_analyses; //!< Activated analyses
    std::map<unsigned, orbital_params> m_oparams;

    orbital_type::flag_t m_ot; //!< Orbital types activated for export
    density_type::flag_t m_dt; //!< Density types activated for export

    enum export_type m_export_dens; //!< How to export densities
    enum export_type m_export_orbs; //!< How to export orbitals

//    std::auto_ptr<molcas_export_h5orb> m_ocore; //!< Pointer to orbital export core
    std::auto_ptr<base_data> m_moldata; //!< Molecular data
    
public:
    /** \brief Constructor
     **/
    molcas_wf_analysis_data(H5::H5File &file);

    /** \brief Virtual destructor
     **/
    virtual ~molcas_wf_analysis_data() { cleanup(); }
    
    /** \brief Initialize orbital export
        \param oe How to export orbitals (possible values: cube, molden)
        \param ot Types of orbitals to export
     **/
    void init_orbital_export(const std::string &oe,
            const orbital_type::flag_t &ot);

    /** \brief Initialize population analysis
        \param name Name of population analysis (possible values: mulliken,
            loewdin)
     **/
    void init_pop_analysis(const std::string &name);

    /** \brief Initialize CT number analysis
        \param name Name of CT number analysis (possible values: atomic)
     **/
    void init_ctnum_analysis(const std::string &name);

    /** \brief Activate certain analysis
     **/
    void activate(enum analysis_type t);
    
    /** \brief Set orbital parameters for one orbital type
        \param t Orbital type
        \param nno Number of leading orbitals
        \param thresh Threshold of important orbitals
     **/
    void set_orbital_params(enum orbital_type::ot t,
            size_t nno, double thresh);
    
    /** \brief Return if a certain analysis should be performed
     **/
    bool is_active(enum analysis_type t);
    
    /** \brief Return parameters for orbital analyses
        \return Pair comprising the number of leading orbitals to print and
            an threshold for important orbitals (see e.g. \ref no_analysis
            for details)
     **/
    orbital_params get_orbital_params(enum orbital_type::ot t);

    /** \brief Retrieve the AO overlap matrix
     **/
    const arma::mat &overlap() {
        return m_moldata->s;
    }

    /** \brief Retrieve the MO coefficient matrices
     **/
    const ab_matrix &coefficients() {
        return m_moldata->c_fb;
    }

    /** \brief Construct a printer of density matrices
        \param name Name of state to which the density matrices belong
            (should be usable as file name)
        \param desc Description of state
        \return Pointer to new density printer
     **/
    density_printer_i *density_printer(const std::string &name,
            const std::string &desc);

    /** \brief Construct a printer of orbitals
        \param name Name of state to which the orbitals belong
            (should be usable as file name)
        \param desc Description of state
        \return Pointer to new orbital printer
     **/
    orbital_printer_i *orbital_printer(const std::string &name,
            const std::string &desc);

    //! \name Population analysis related functions
    //@{

    /** \brief Number of population analyses available
     **/
    size_t n_pop_analyses() { return m_pa.size(); }

    /** \brief Name of i-th population analysis
     **/
    const std::string &pop_name(size_t i) {
        return m_pa[i]->name;
    }

    /** \brief Row labels for i-th population analysis
     **/
    const std::vector<std::string> &pop_labels(size_t i) {
        return m_pa[i]->labels;}

    /** \brief i-th population analysis
     **/
    const pop_analysis_i &pop_analysis(size_t i) {
        return *m_pa[i]->analysis;
    }

    /** \brief Base populations for i-th population analysis
     **/
    const arma::vec &ref_population(size_t i) {
        return m_pa[i]->p0;

    }

    //@}

    //! \name CT number analysis related functions
    //@{

    /** \brief Number of CT number analyses available
     **/
    size_t n_ctnum_analyses() {}

    /** \brief Name of i-th CT number analysis
     **/
    const std::string &ctnum_name(size_t i) {}

    /** \brief i-th CT number analysis
     **/
    const ctnum_analysis_i &ctnum_analysis(size_t i) {}

    /** \brief Printer of i-th CT number data
     **/
    std::auto_ptr<ctnum_printer_i> ctnum_printer(size_t i,
            const std::string &name, const std::string &desc) {}

    //@}

    /** \brief Builder of exciton moments
     **/
    const mom_builder_i &mom_builder() {
        return m_moldata->mom;
    }

    /** \brief Build the density matrix
        \param buf Buffer with density matrix data
        \return Full density matrix in the AO basis
        
        The appropriate number of doubly occupied orbitals are added
        and the density is transformed to the AO basis.
     **/
    ab_matrix build_dm(const double *buf, const double *sbuf, const bool aeqb_dens);
    
private:
    void initialize();
    void cleanup();
    std::string get_mo_types(const H5::DataSet &Set);
    arma::mat get_mo_vectors(const H5::DataSet &Set, const size_t nsym, const int *nbas);
};

/** \brief Setup the wave function / density matrix analysis data for Molcas
    \return Pointer to data object

    This is a convenience wrapper to initialize the analysis data object. If
    more fine-tuning is required, the function can serve as template on how
    the object can be setup.

    \ingroup libwfa_molcas
 **/
molcas_wf_analysis_data *molcas_setup_wf_analysis_data(H5::H5File file);

} // libwfa

#endif // LIBWFA_MOLCAS_WF_ANALYSIS_DATA_H