#ifndef LIBWFA_MOLCAS_WF_ANALYSIS_DATA_H
#define LIBWFA_MOLCAS_WF_ANALYSIS_DATA_H

#include <libwfa/wf_analysis_data_i.h>
namespace libwfa
{
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
//        molcas_mom_builder mom; //!< Moments builder

        /** \brief Constructor
            \param nao Number of atomic basis functions
            \param nb2 Number of shell pairs
            \param maxmm Maximum multipole moment available
            \param unr Unrestricted calculation
         **/
        base_data(size_t nao, size_t nb2, size_t maxmm, bool unr) :
            c_fb(unr) { } // add mom ...
    };

    /** \brief Export types
     **/
    enum export_type {
        EXPORT_NONE = 0,//!< EXPORT_NONE
        EXPORT_MOLDEN,  //!< EXPORT_MOLDEN
    };
    
private:
    std::vector<pa_data *> m_pa; //!< Population analyses
//    std::vector<cta_data *> m_cta; //!< CT number analyses

    std::bitset<WFA_TYPES> m_analyses; //!< Activated analyses
    std::map<unsigned, orbital_params> m_oparams;

    orbital_type::flag_t m_ot; //!< Orbital types activated for export
    density_type::flag_t m_dt; //!< Density types activated for export

//    enum export_type m_export_dens; //!< How to export densities
//    enum export_type m_export_orbs; //!< How to export orbitals

//    std::auto_ptr<qchem_export_cube> m_ccore; //!< Pointer to cube export core
//    std::auto_ptr<qchem_export_molden> m_mcore; //!< Pointer to molden export core
    std::auto_ptr<base_data> m_moldata; //!< Molecular data
    
public:
    /** \brief Constructor
     **/
    molcas_wf_analysis_data();

    /** \brief Virtual destructor
     **/
    virtual ~molcas_wf_analysis_data();    
};
} // libwfa
    
#endif // LIBWFA_MOLCAS_WF_ANALYSIS_DATA_H