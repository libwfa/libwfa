#include <libwfa/libwfa.h>
#include "molcas_wf_analysis_data.h"

namespace libwfa {
using namespace H5;

const char molcas_wf_analysis_data::k_clazz[] = "molcas_wf_analysis_data";
    
molcas_wf_analysis_data::molcas_wf_analysis_data(H5File &file) :
    m_file(file), m_export_dens(EXPORT_NONE), m_export_orbs(EXPORT_NONE) {
    initialize();
    std::cout << "Initialization finished." << std::endl;
}

void molcas_wf_analysis_data::init_pop_analysis(const std::string &name) {

    const std::vector<std::string> &atoms = m_moldata->atoms;
    const arma::mat &s = m_moldata->s;
    const arma::uvec &b2a = m_moldata->bf2atoms;
    const arma::vec &p0 = m_moldata->atomic_charges;

    if (name == "mulliken") {
        m_pa.push_back(new pa_data("Mulliken Population Analysis",
                atoms, new pop_mulliken(s, b2a), p0));
    }
    else if (name == "loewdin") {
        m_pa.push_back(new pa_data("Loewdin Population Analysis",
                atoms, new pop_loewdin(s, b2a), p0));
    }
}

void molcas_wf_analysis_data::init_ctnum_analysis(const std::string &name) {

    /*
    //const arma::mat &s = m_moldata->mom.get(0, 0);
    arma::mat s;
    const arma::uvec &b2a = m_moldata->bf2atoms;
    /*if (name  == "atomic") {
        m_cta.push_back(new cta_data("Atomic CT numbers", "atomic",
                new libwfa::ctnum_analysis(s, b2a)));
    }*/
}

void molcas_wf_analysis_data::activate(enum analysis_type t) {

    if (t < 0 || t > (int)WFA_TYPES) return;
    m_analyses.set(t, true);
}

void molcas_wf_analysis_data::set_orbital_params(enum orbital_type::ot t,
        size_t nno, double thresh) {

    if ((t & (orbital_type::NO | orbital_type::NDO | orbital_type::NTO)) == 0)
        return;

    m_oparams[t] = orbital_params(nno, thresh);
}

molcas_wf_analysis_data::orbital_params
molcas_wf_analysis_data::get_orbital_params(enum orbital_type::ot t) {

    std::map<unsigned, orbital_params>::const_iterator i = m_oparams.find(t);
    if (i == m_oparams.end()) return orbital_params();
    else return i->second;
}

bool molcas_wf_analysis_data::is_active(enum analysis_type t) {

    if (t < 0 || t > (int)WFA_TYPES) return false;
    return m_analyses[t];
}

density_printer_i *molcas_wf_analysis_data::density_printer(
    const std::string &name, const std::string &desc) {

    /*if (m_export_dens == EXPORT_FCHK)
        return new qchem_density_printer_fchk(desc, m_dt);
    else if (m_export_dens == EXPORT_CUBE)
        return new density_printer_cube(*m_ccore, name, desc, m_dt);
    else*/
        return new density_printer_nil();
}

orbital_printer_i *molcas_wf_analysis_data::orbital_printer(
        const std::string &name, const std::string &desc) {

/*    if (m_export_orbs == EXPORT_MOLDEN)
        return new orbital_printer_molden(*m_mcore, name, m_ot);
    else if (m_export_orbs == EXPORT_CUBE)
        return new orbital_printer_cube(*m_ccore, name, desc, m_ot);
    else */
        return new orbital_printer_nil();
}

ab_matrix molcas_wf_analysis_data::build_dm(const double *buf) {
    int nao = m_moldata->c_fb.nrows_a();
    int nmo = m_moldata->c_fb.ncols_a();
    
    ab_matrix dmo = ab_matrix(nmo, nmo);
    dmo.alpha().zeros();
    
    int nocc = m_moldata->ndocc;
    for (int iocc=0; iocc<nocc; iocc++)
        dmo.alpha().at(iocc, iocc) = 1.;
        
    int nact = m_moldata->nact;
    arma::mat dmact(buf, nact, nact);
    dmo.alpha().submat(nocc, nocc, nocc+nact-1, nocc+nact-1) = 0.5 * dmact;

    ab_matrix dao(nao, nao);
    dao.alpha() = m_moldata->c_fb.alpha() * dmo.alpha() * m_moldata->c_fb.alpha().t();
        
    return dao;
}

void molcas_wf_analysis_data::initialize() {
    
    static const char method[] = "initialize()";    
    
    Group Grp_main = m_file.openGroup("/");
    hsize_t natoms;
    int nbas_t;
    hsize_t nmo;
    
    // Read basic info
    {
        Attribute Att_nbas = Grp_main.openAttribute("NBAS");
        DataSpace Spac_nbas = Att_nbas.getSpace();
        hsize_t nnbas;
        Spac_nbas.getSimpleExtentDims(&nnbas, NULL);
        
        int nbas[nnbas];
        Att_nbas.read(PredType::NATIVE_INT, nbas);
    
        nbas_t = 0;
        std::cout << " Basis functions per irrep: ";
        for (int ibas=0; ibas<nnbas; ibas++) {
            std::cout << " " << nbas[ibas];
            nbas_t += nbas[ibas];
        }
        std::cout << std::endl;    
    
        m_moldata = std::auto_ptr<base_data>(new base_data(nbas_t, 0, 0, false));
    }
    
    // Read atomic numbers/charges. Different in the case of ECPs(?)
    {
        DataSet Set = m_file.openDataSet("CENTER_CHARGES");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 1)         
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent rank");
        
        Space.getSimpleExtentDims(&natoms, NULL);
        std::cout << " Number of atoms: " << natoms << std::endl;
        //double t[natoms];
        int num_buf[natoms];
        Set.read(&num_buf, PredType::NATIVE_INT);
        
        m_moldata->atoms = std::vector<std::string>(natoms);
        m_moldata->atomic_numbers = arma::uvec(natoms);
        m_moldata->atomic_charges = arma::vec(natoms);
        for (size_t i = 0; i < natoms; i++) {
            m_moldata->atomic_numbers.at(i) = num_buf[i];
            m_moldata->atomic_charges.at(i) = num_buf[i];
        }
        std::cout << std::endl;
    }
    
    // Read atom labels
    {
        DataSet Set = m_file.openDataSet("CENTER_LABELS");
        int len = 6;
        StrType strtype = Set.getStrType();      
        char labels[natoms][len];
        Set.read(labels, strtype);
        
        for (size_t i = 0; i < natoms; i++) {
            m_moldata->atoms.at(i) = std::string(labels[i], len);
        }
    }
    
    // Read atomic coordinates
    {
        DataSet Set = m_file.openDataSet("CENTER_COORDINATES");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 2) {            
            throw libwfa_exception(k_clazz,
                method, __FILE__, __LINE__, "Inconsistent rank");
        }
        
        hsize_t dims[rank];
        Space.getSimpleExtentDims(dims, NULL);
        double buf[dims[0] * dims[1]];
        Set.read(&buf, PredType::NATIVE_DOUBLE);
        m_moldata->coordinates = arma::mat(buf, 3, natoms);
        //m_moldata->coordinates.print();
    }
    
    // MO types
    {
        DataSet Set = m_file.openDataSet("MO_TYPEINDICES");
        int len = 1;
        StrType strtype = Set.getStrType();      

        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 1) 
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent rank");
        
        Space.getSimpleExtentDims(&nmo, NULL);
        char buf[nmo][len];
        Set.read(buf, strtype);
        
        std::string str(*buf, nmo*len);
        m_moldata->ndocc = std::count(str.begin(), str.end(), 'I');
        m_moldata->nact  = std::count(str.begin(), str.end(), '1') +
                           std::count(str.begin(), str.end(), '2') +
                           std::count(str.begin(), str.end(), '3');
    }
    
    // Mapping of basis functions to atoms
    //   Use DESYMM in the case of symmetry
    {
        DataSet Set = m_file.openDataSet("BASIS_FUNCTION_IDS");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 2) 
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent rank");
        
        hsize_t dims[rank];
        Space.getSimpleExtentDims(dims, NULL);
        if (dims[1] != 4)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent dims");
        
        int buf[dims[0] * dims[1]];
        Set.read(&buf, PredType::NATIVE_INT);
        m_moldata->bf2atoms = arma::uvec(dims[0]);
        for (size_t i = 0; i < dims[0]; i++) {
            m_moldata->bf2atoms.at(i) = buf[4*i] - 1;
        }
        //m_moldata->bf2atoms.print();
    }
    
    // MO-coefficients
    //  TODO: transpose?
    {
        DataSet Set = m_file.openDataSet("MO_VECTORS");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 1) 
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent rank");
        
        hsize_t dim;
        Space.getSimpleExtentDims(&dim, NULL);
        if (nbas_t * nbas_t != dim)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent MO vectors");
        
        double buf[dim];
        Set.read(&buf, PredType::NATIVE_DOUBLE);
        m_moldata->c_fb.alpha() = arma::mat(buf, nbas_t, nbas_t);
    }
    
    // Overlap matrix
    {
        DataSet Set = m_file.openDataSet("AO_OVERLAP_MATRIX");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 1) 
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent rank");
        
        hsize_t dim;
        Space.getSimpleExtentDims(&dim, NULL);
        if (nbas_t * nbas_t != dim)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent S matrix");
        
        double buf[dim];
        Set.read(&buf, PredType::NATIVE_DOUBLE);
        m_moldata->s = arma::mat(buf, nbas_t, nbas_t);
    }
}

void molcas_wf_analysis_data::cleanup() {

    delete m_moldata.release();

    for (std::vector<pa_data *>::iterator i = m_pa.begin();
            i != m_pa.end(); i++) {
        delete *i; *i = 0;
    }
    m_pa.clear();
/*    for (std::vector<cta_data *>::iterator i = m_cta.begin();
            i != m_cta.end(); i++) {
        delete *i; *i = 0;
    }
    m_cta.clear();*/
}

molcas_wf_analysis_data *molcas_setup_wf_analysis_data(H5::H5File file) {
    std::cout << "Starting setup..." << std::endl;
    
    molcas_wf_analysis_data *h = new molcas_wf_analysis_data(file);
    
    size_t norb = 3;
    double thresh = 1e-3;
    
    // Setup parameters for orbital print-out; currently use the same parameters
    // for NOs, NDOs, and NTOs
    h->set_orbital_params(orbital_type::NO,  norb, thresh);
    h->set_orbital_params(orbital_type::NDO, norb, thresh);
    h->set_orbital_params(orbital_type::NTO, norb, thresh);    
    
    // Activate Mulliken population analysis
    // TODO: Loewdin
    h->init_pop_analysis("mulliken");
    
    h->activate(molcas_wf_analysis_data::NO);
    h->activate(molcas_wf_analysis_data::NDO);
    h->activate(molcas_wf_analysis_data::FORM_AD);
    
    return h;
}
} // namespace libwfa