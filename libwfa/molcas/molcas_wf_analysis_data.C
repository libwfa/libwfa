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

ab_matrix molcas_wf_analysis_data::build_dm(const double *buf, const double *sbuf, const bool aeqb_dens) {
    int nao = m_moldata->c_fb.nrows_a();
    int nmo = m_moldata->c_fb.ncols_a();

    bool aeqb = m_moldata->c_fb.is_alpha_eq_beta() && aeqb_dens;

    ab_matrix dmo(aeqb);
    { // Construct the alpha part
        dmo.alpha() = arma::zeros(nmo, nmo);
        if (!aeqb) dmo.beta() = arma::zeros(nmo, nmo);
        std::string imot, jmot, imot_b;
        for (int imo=0; imo<nmo; imo++) {
            imot = m_moldata->mo_types_a[imo];
            imot_b = m_moldata->mo_types_b[imo];
            if (imot == "I" || imot == "F"){ // Inactive or frozen
                dmo.alpha().at(imo, imo) = 1.;
            }
            else if (imot == "1" || imot == "2" || imot == "3") { // Active
                if (imot != imot_b)
                    throw libwfa_exception(k_clazz, "build_dm", __FILE__, __LINE__, "Inconsistent mo_types");
                
                for (int jmo=0; jmo<nmo; jmo++){
                    jmot = m_moldata->mo_types_a[jmo];
                    if (jmot == "1" || jmot == "2" || jmot == "3") {
                        dmo.alpha().at(imo, jmo) = *buf * 0.5;
                        
                        if (!aeqb_dens) {
                            dmo.alpha().at(imo, jmo) += *sbuf * 0.5;
                            dmo.beta().at(imo, jmo)  =(*buf - *sbuf)*0.5;
                        }
                        buf++; sbuf++;
                    }
                }
            }
            else if (imot == "S") {} // Secondary
            else {
                std::ostringstream os;
                os << std::endl << "Unknown MO type: " << imot;
                const std::string errmsg = os.str();
                throw libwfa_exception(k_clazz, "build_dm", __FILE__, __LINE__, errmsg.c_str());
            }
            
            if (!aeqb) {
                if (imot_b == "I" || imot_b == "F")
                    dmo.beta().at(imo, imo) = 1.;
            }
        }
    }

    ab_matrix dao(aeqb);
    dao.alpha() = m_moldata->c_fb.alpha() * dmo.alpha() * m_moldata->c_fb.alpha().t();
    if (!aeqb)
        dao.beta() = m_moldata->c_fb.beta() * dmo.beta() * m_moldata->c_fb.beta().t();

    return dao;
}

void molcas_wf_analysis_data::initialize() {

    static const char method[] = "initialize()";

    hsize_t natoms;
    hsize_t nsym, nbas_t;
    int nbas[8];
    arma::mat desym;
    bool aeqb;

    Group Grp_main = m_file.openGroup("/");

    try {
        DataSet Set = m_file.openDataSet("MO_TYPEINDICES");
        aeqb = true;
    }
    catch( FileIException error ) { // Unrestricted case
        std::cout << std::endl << "MO_TYPEINDICES not found" << std::endl;
        std::cout << "  ... assuming unrestricted orbitals." << std::endl << std::endl;
        aeqb = false;
    }

    // Number of basis functions
    {
        Attribute Att_nbas = Grp_main.openAttribute("NBAS");
        DataSpace Spac_nbas = Att_nbas.getSpace();
        Spac_nbas.getSimpleExtentDims(&nsym, NULL);
        if (nsym > 8)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "More than 8 irreps");
        //m_moldata->nsym = nsym;

        Att_nbas.read(PredType::NATIVE_INT, nbas);

        nbas_t = 0;
        std::cout << " Basis functions per irrep: ";
        for (size_t isym=0; isym<nsym; isym++) {
            std::cout << " " << nbas[isym];
            nbas_t += nbas[isym];
        }
        std::cout << std::endl;

        m_moldata = std::auto_ptr<base_data>(new base_data(nbas_t, 0, 0, aeqb));
    }

    // Read atomic numbers/charges. Different in the case of ECPs(?)
    {
        H5std_string h5label = (nsym==1) ? "CENTER_CHARGES" : "DESYM_CENTER_CHARGES";
        DataSet Set = m_file.openDataSet(h5label);
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 1)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent rank");

        Space.getSimpleExtentDims(&natoms, NULL);
        std::cout << " Number of atoms: " << natoms << std::endl;
        int num_buf[natoms];
        Set.read(&num_buf, PredType::NATIVE_INT);

        m_moldata->atoms = std::vector<std::string>(natoms);
        m_moldata->atomic_numbers = arma::uvec(natoms);
        m_moldata->atomic_charges = arma::vec(natoms);
        for (size_t i = 0; i < natoms; i++) {
            m_moldata->atomic_numbers.at(i) = num_buf[i];
            m_moldata->atomic_charges.at(i) = num_buf[i];
        }
    }

    // Read atom labels
    {
        H5std_string h5label;
        hsize_t len;
        if (nsym==1) {
            h5label = "CENTER_LABELS";
            len = 6;
        }
        else {
            h5label = "DESYM_CENTER_LABELS";
            len = 10;
        }
        DataSet Set = m_file.openDataSet(h5label);

        char labels[natoms][len];

        StrType strtype = Set.getStrType();
        Set.read(labels, strtype);

        for (size_t i = 0; i < natoms; i++) {
            m_moldata->atoms.at(i) = std::string(labels[i], len);
        }
    }

    // Read atomic coordinates
    {
        H5std_string h5label = (nsym==1) ? "CENTER_COORDINATES" : "DESYM_CENTER_COORDINATES";
        DataSet Set = m_file.openDataSet(h5label);
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
        m_moldata->coordinates.print();
    }

    // Desymmetrization matrix
    if (nsym>1){
        H5std_string h5label = "DESYM_MATRIX";
        DataSet Set = m_file.openDataSet(h5label);
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 1)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent rank");

        hsize_t dim;
        Space.getSimpleExtentDims(&dim, NULL);
        double buf[dim];
        Set.read(&buf, PredType::NATIVE_DOUBLE);

        if (nbas_t * nbas_t != dim)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent desym matrix");

        desym = arma::mat(buf, nbas_t, nbas_t);
    }

    // Mapping of basis functions to atoms
    {
        H5std_string h5label = (nsym==1) ? "BASIS_FUNCTION_IDS" : "DESYM_BASIS_FUNCTION_IDS";
        DataSet Set = m_file.openDataSet(h5label);
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
        double buf[dim];
        Set.read(&buf, PredType::NATIVE_DOUBLE);


        if (nsym == 1) { // no symmetry
            if (nbas_t * nbas_t != dim)
                throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent S-matrix (no symmetry)");

                m_moldata->s = arma::mat(buf, nbas_t, nbas_t);
        } // symmetry
        else {
            hsize_t dim_chk = 0;
            for (size_t isym=0; isym<nsym; isym++)
                dim_chk += nbas[isym] * nbas[isym];

            if (dim_chk != dim)
                throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent S-matrix (symmetry)");

            // Fill the blocks
            m_moldata->s = arma::zeros(nbas_t, nbas_t);
            size_t i=0;
            double *buf_ptr = buf;
            for (size_t isym=0; isym<nsym; isym++) {
                int kbas = nbas[isym];
                if (kbas==0) continue;

                m_moldata->s.submat(i, i, i+kbas-1, i+kbas-1) = arma::mat(buf_ptr, kbas, kbas);

                i+=kbas;
                buf_ptr += kbas*kbas;
            }
//            std::cout << "fptmp: s" << std::endl;
//            m_moldata->s.print();
            m_moldata->s = desym * m_moldata->s * desym.t();
        }
    }

    // MO types
    if (aeqb) {
        DataSet Set = m_file.openDataSet("MO_TYPEINDICES");
        std::string str = get_mo_types(Set);
        m_moldata->mo_types_a = str;
        m_moldata->mo_types_b = str;
    }
    else { // Unrestricted case
        {
            DataSet Set = m_file.openDataSet("MO_ALPHA_TYPEINDICES");
            m_moldata->mo_types_a = get_mo_types(Set);
        }
        {
            DataSet Set = m_file.openDataSet("MO_BETA_TYPEINDICES");
            m_moldata->mo_types_b = get_mo_types(Set);
        }
        std::cout << m_moldata->mo_types_b << std::endl;
    }

    // MO-coefficients
    // TODO: adjust dimensions in case of rectangular MO-matrix
    if (aeqb) {
        DataSet Set = m_file.openDataSet("MO_VECTORS");
        m_moldata->c_fb.alpha() = get_mo_vectors(Set, nsym, nbas);
        if (nsym>1)
            m_moldata->c_fb.alpha() = desym * m_moldata->c_fb.alpha();
    }
    else {
        {
            DataSet Set = m_file.openDataSet("MO_ALPHA_VECTORS");
            m_moldata->c_fb.alpha() = get_mo_vectors(Set, nsym, nbas);
            if (nsym>1)
                m_moldata->c_fb.alpha() = desym * m_moldata->c_fb.alpha();
        }
        {
            DataSet Set = m_file.openDataSet("MO_BETA_VECTORS");
            m_moldata->c_fb.beta() = get_mo_vectors(Set, nsym, nbas);
            if (nsym>1)
                m_moldata->c_fb.beta() = desym * m_moldata->c_fb.beta();
        }
    }
}

std::string molcas_wf_analysis_data::get_mo_types(const H5::DataSet &Set) {
    int len = 1;
    StrType strtype = Set.getStrType();

    DataSpace Space = Set.getSpace();
    int rank = Space.getSimpleExtentNdims();
    if (rank != 1)
        throw libwfa_exception(k_clazz, "get_mo_types", __FILE__, __LINE__, "Inconsistent rank");

    hsize_t nmo;
    Space.getSimpleExtentDims(&nmo, NULL);
    char buf[nmo][len];
    Set.read(buf, strtype);

    std::string str(*buf, nmo*len);
    return str;
}

arma::mat molcas_wf_analysis_data::get_mo_vectors(const H5::DataSet &Set, const size_t nsym, const int *nbas) {
    DataSpace Space = Set.getSpace();
    int rank = Space.getSimpleExtentNdims();
    if (rank != 1)
        throw libwfa_exception(k_clazz, "get_mo_vectors", __FILE__, __LINE__, "Inconsistent rank");

    hsize_t dim;
    Space.getSimpleExtentDims(&dim, NULL);
    double buf[dim];
    Set.read(&buf, PredType::NATIVE_DOUBLE);

    int nbas_t = 0;
    for (size_t isym=0; isym<nsym; isym++)
        nbas_t += nbas[isym];

    arma::mat retmat;
    if (nsym == 1) { // no symmetry
        if (nbas_t * nbas_t != dim)
            throw libwfa_exception(k_clazz, "get_mo_vectors", __FILE__, __LINE__, "Inconsistent MO vectors (no symmetry)");

        retmat = arma::mat(buf, nbas_t, nbas_t);
    } // symmetry
    else {
        hsize_t dim_chk = 0;
        for (size_t isym=0; isym<nsym; isym++)
            dim_chk += nbas[isym] * nbas[isym];

        if (dim_chk != dim)
            throw libwfa_exception(k_clazz, "get_mo_vectors", __FILE__, __LINE__, "Inconsistent MO vectors (symmetry)");

        // Fill the blocks of the MO matrix
        retmat = arma::zeros(nbas_t, nbas_t);
        size_t i=0;
        double *buf_ptr = buf;
        for (size_t isym=0; isym<nsym; isym++) {
            int kbas = nbas[isym];
            if (kbas==0) continue;

            retmat.submat(i, i, i+kbas-1, i+kbas-1) = arma::mat(buf_ptr, kbas, kbas);

            i+=kbas;
            buf_ptr += kbas*kbas;
        }
    }
    return retmat;
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