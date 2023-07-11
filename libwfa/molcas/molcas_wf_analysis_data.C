//************************************************************************
//* This file is part of libwfa.                                         *
//*                                                                      *
//* libwfa is free software; you can redistribute and/or modify          *
//* it under the terms of the BSD 3-Clause license.                      *
//* libwfa is distributed in the hope that it will be useful, but it     *
//* is provided "as is" and without any express or implied warranties.   *
//* For more details see the full text of the license in the file        *
//* LICENSE.                                                             *
//*                                                                      *
//* Copyright (c) 2014, F. Plasser and M. Wormit. All rights reserved.   *
//* Modifications copyright (C) 2019, Loughborough University.           *
//************************************************************************


#include <fstream>
#include <string>
#include <libwfa/libwfa.h>
#include <libwfa/core/constants.h>
#include "molcas_wf_analysis_data.h"

namespace libwfa {
using namespace H5;

const char molcas_wf_analysis_data::k_clazz[] = "molcas_wf_analysis_data";

molcas_wf_analysis_data::molcas_wf_analysis_data(char *inp) :
    m_export_dens(EXPORT_NONE), m_export_orbs(EXPORT_NONE) {
    read_input(inp);

    std::cout << "Opening Molcas HDF5 file " << m_input->file_name << std::endl;
    m_file = H5File(m_input->file_name, H5F_ACC_RDWR);

    initialize();
    std::cout << "Initialization finished." << std::endl;
}

void molcas_wf_analysis_data::init_orbital_export(const std::string &oe,
        const orbital_type::flag_t &ot) {

    m_ot = ot;
    if (oe == "h5" || oe == "hdf5") { m_export_orbs = EXPORT_H5; setup_h5core(); }
    else { m_ot.reset(); }
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
    else if (name == "lowdin") {
        m_pa.push_back(new pa_data("Lowdin Population Analysis",
                atoms, new pop_loewdin(s, b2a), p0));
    }
    else {
        throw libwfa_exception(k_clazz, "build_dm", __FILE__, __LINE__, "Unknown population analysis type.");
    }
}

void molcas_wf_analysis_data::init_ctnum_analysis(const std::string &name) {

    const arma::mat &s = m_moldata->s;
    const arma::uvec &b2a = m_moldata->bf2atoms;

    const std::vector<std::string> &prop_list = m_input->prop_list;
    const std::vector<std::vector<int>> &at_lists = m_input->at_lists;


    if (name  == "mulliken") {
        m_cta.push_back(new cta_data("CT numbers (Mulliken)", "mulliken",
                new libwfa::ctnum_analysis(s, b2a, name, prop_list, at_lists)));
    }
    else if (name == "lowdin") {
        m_cta.push_back(new cta_data("CT numbers (Lowdin)", "lowdin",
                new libwfa::ctnum_analysis(s, b2a, name, prop_list, at_lists)));
    }
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

    return new density_printer_nil();
}

orbital_printer_i *molcas_wf_analysis_data::orbital_printer(
        const std::string &name, const std::string &desc) {
    if (m_export_orbs == EXPORT_H5)
        return new orbital_printer_molden(*m_h5core, name, m_ot);
    else
        return new orbital_printer_nil();
}

std::unique_ptr<ctnum_printer_i> molcas_wf_analysis_data::ctnum_printer(size_t i,
            const std::string &name, const std::string &desc) {
    ctnum_printer_i *pr = new ctnum_export(name + "_ctnum_" +
            m_cta[i]->suffix, desc);
    return std::unique_ptr <ctnum_printer_i>(pr);
}

ab_matrix molcas_wf_analysis_data::build_dm(const double *buf, const double *sbuf, const bool aeqb_dens) {
    //int nao = m_moldata->c_fb.nrows_a();
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

ab_matrix molcas_wf_analysis_data::build_dm_ao(const double *buf, const double *sbuf, const size_t dim) {
    // This is the default, a totally symmetric density matrix.
    return build_dm_ao(buf, sbuf, dim, 1, 1);
}

ab_matrix molcas_wf_analysis_data::build_dm_ao(const double *buf, const double *sbuf, const size_t dim,
    const int isym, const int jsym) {

    int nao = m_moldata->c_fb.nrows_a();
    size_t nsym = m_moldata->nbas.size();

    int psym = m_moldata->symmult(isym-1, jsym-1);

    arma::mat den, sden;
    read_ao_mat(buf, dim, den, nsym, psym);

    if (sbuf == NULL)
        sden = arma::zeros(nao, nao);
    else
        read_ao_mat(sbuf, dim, sden, nsym, psym);

    bool aeqb = true;
    {
        double dnorms = accu(den%den);
        double sdnorms = accu(sden%sden);

        if ((dnorms > 1.e-8) && (sdnorms > 1.e-8)) {
            aeqb = false;
        }
    }

    ab_matrix dao(aeqb);
    dao.alpha() = 0.5 * (den + sden);
    if (!aeqb)
        dao.beta() = 0.5 * (den - sden);

    if (nsym > 1) {
        dao.alpha() = m_moldata->desym * dao.alpha() * m_moldata->desym.t();
        if (!aeqb)
            dao.beta() = m_moldata->desym * dao.beta() * m_moldata->desym.t();
    }

    return dao;
}

arma::vec molcas_wf_analysis_data::read_vec_h5(H5std_string key) {
    DataSet Set = m_file.openDataSet(key);
    DataSpace Space = Set.getSpace();
    if (Space.getSimpleExtentNdims() != 1)
        throw libwfa_exception(k_clazz, "read_vec_h5", __FILE__, __LINE__, "Inconsistent rank");

    hsize_t dim;
    Space.getSimpleExtentDims(&dim, NULL);

    arma::vec retvec = arma::vec(dim);
    double *buf = retvec.memptr();
    Set.read(buf, PredType::NATIVE_DOUBLE);

    return retvec;
}

arma::cube molcas_wf_analysis_data::read_cube_h5(H5std_string key) {
    DataSet Set = m_file.openDataSet(key);
    DataSpace Space = Set.getSpace();
    if (Space.getSimpleExtentNdims() != 3)
        throw libwfa_exception(k_clazz, "read_cube_h5", __FILE__, __LINE__, "Inconsistent rank");

    hsize_t dims[3];
    Space.getSimpleExtentDims(dims, NULL);

    arma::cube dens = arma::cube(dims[2], dims[1], dims[0]);
    double *dens_buf = dens.memptr();
    Set.read(dens_buf, PredType::NATIVE_DOUBLE);

    return dens;
}

void molcas_wf_analysis_data::energy_print(const double ener, std::ostream &out) {
    out << "Energy: " << std::setprecision(5) << ener*constants::au2eV
        << " eV" << std::endl << std::endl;
}

std::string molcas_wf_analysis_data::rasscf_label() {

    Group Grp_main = m_file.openGroup("/");

    int lsym;
    {
        Attribute Att = Grp_main.openAttribute("LSYM");
        Att.read(PredType::NATIVE_INT, &lsym);
    }
    int imult;
    {
        Attribute Att = Grp_main.openAttribute("SPINMULT");
        Att.read(PredType::NATIVE_INT, &imult);
    }

    std::ostringstream ostream;

    std::string str = m_moldata->irrep_labels.at(lsym-1);
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);

    ostream << "(" << imult << ") ";
    ostream << str;

    return ostream.str();
}

std::vector<std::string> molcas_wf_analysis_data::rassi_labels(int *mult, int *irrep, int nstate) {
    std::vector<std::string> state_labels;

    int offset;
    char lab [] = {'S', 'D', 'T', 'Q', 'Q', 'H', 'H'};
    int  ind [] = { 0 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 };

    Group Grp_main = m_file.openGroup("/");

    {
        Attribute Att = Grp_main.openAttribute("STATE_SPINMULT");
        Att.read(PredType::NATIVE_INT, mult);
    }
    {
        Attribute Att = Grp_main.openAttribute("STATE_IRREPS");
        Att.read(PredType::NATIVE_INT, irrep);
    }

    for (int istate = 0; istate < nstate; istate++) {
        offset = mult[istate] - 1;
        std::ostringstream name;
        name << lab[offset] << ind[offset]++;
        state_labels.push_back(name.str());
    }

    return state_labels;
}

void molcas_wf_analysis_data::initialize() {

    static const char method[] = "initialize()";

    hsize_t natoms;
    hsize_t nsym, nbas_t;
    int nbas[8];

    Group Grp_main = m_file.openGroup("/");

    bool aeqb = true;
    bool found_mos = true;

    if (! m_input->debug) Exception::dontPrint(); // Do not print excessive error messages from HDF5

    try {
        DataSet Set = m_file.openDataSet("MO_VECTORS");
        std::cout << std::endl << "Found restricted MO-coefficients: MO_VECTORS" << std::endl;
    }
    catch( FileIException error ) {
        try {
            DataSet Set = m_file.openDataSet("MO_ALPHA_VECTORS");
            aeqb = false;
            std::cout << std::endl << "Found unrestricted MO-coefficients: MO_ALPHA_VECTORS" << std::endl;
        }
        catch( FileIException error ) {
            found_mos = false;
            std::cout << std::endl << "Did not find any MO-coefficients" << std::endl;
        }
    }
    // Exception::setAutoPrint(&H5Eprint, &std::cerr); // How does this work?

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
        for (size_t isym=0; isym<nsym; isym++) {
            nbas_t += nbas[isym];
        }
        m_moldata = std::unique_ptr<base_data>(new base_data(nbas_t, 2, aeqb));

        m_moldata->nbas = arma::uvec(nsym);
        for (size_t isym=0; isym<nsym; isym++) {
            m_moldata->nbas(isym) = nbas[isym];
        }
    }

    // Irrep labels
    {
        Attribute Att = Grp_main.openAttribute("IRREP_LABELS");
        DataSpace Space = Att.getSpace();
        size_t len = 3;

        char labels[nsym][len];

        StrType strtype = Att.getStrType();
        Att.read(strtype, labels);
        m_moldata->irrep_labels = std::vector<std::string>(nsym);
        for (size_t isym = 0; isym < nsym; isym++) {
            m_moldata->irrep_labels.at(isym) = std::string(labels[isym], len);
        }
    }

    // Molcas module
    {
        Attribute Att = Grp_main.openAttribute("MOLCAS_MODULE");
        StrType strtype = Att.getStrType();
        Att.read(strtype, m_moldata->molcas_module);
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
        //m_moldata->coordinates.print();
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
        double* buf = new double[dim];
        Set.read(buf, PredType::NATIVE_DOUBLE);

        if (nbas_t * nbas_t != dim)
            throw libwfa_exception(k_clazz, method, __FILE__, __LINE__, "Inconsistent desym matrix");

        m_moldata->desym = arma::mat(buf, nbas_t, nbas_t);
    }

    // Product table of the irreps. This is hardcoded following the pegamoid code.
    {
        m_moldata->symmult = {
          {0, 1, 2, 3, 4, 5, 6, 7},
          {1, 0, 3, 2, 5, 4, 7, 6},
          {2, 3, 0, 1, 6, 7, 4, 5},
          {3, 2, 1, 0, 7, 6, 5, 4},
          {4, 5, 6, 7, 0, 1, 2, 3},
          {5, 4, 7, 6, 1, 0, 3, 2},
          {6, 7, 4, 5, 2, 3, 0, 1},
          {7, 6, 5, 4, 3, 2, 1, 0},
      };
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
        double* buf = new double[dim];
        Set.read(buf, PredType::NATIVE_DOUBLE);

        read_ao_mat(buf, dim, m_moldata->s, nsym, 0);
        if (nsym > 1)
            m_moldata->s = m_moldata->desym * m_moldata->s * m_moldata->desym.t();
    }

    // MO-coefficients
    if (found_mos) {
        // MO types
        if (aeqb) {
            std::string str = get_mo_types("MO_TYPEINDICES");
            m_moldata->mo_types_a = str;
            m_moldata->mo_types_b = str;
        }
        else { // Unrestricted case
            m_moldata->mo_types_a = get_mo_types("MO_ALPHA_TYPEINDICES");
            m_moldata->mo_types_b = get_mo_types("MO_BETA_TYPEINDICES");
        }

        // MO-coefficients
        // TODO: adjust dimensions in case of rectangular MO-matrix
        if (aeqb) {
            m_moldata->c_fb.alpha() = get_mo_vectors("MO_VECTORS");
            //std::cout << "c_fb.alpha" << std::endl;
            //m_moldata->c_fb.alpha().print();
        }
        else {
            m_moldata->c_fb.alpha() = get_mo_vectors("MO_ALPHA_VECTORS");
            m_moldata->c_fb.beta()  = get_mo_vectors("MO_BETA_VECTORS");
        }
    }
    else {
        // Do a Lowdin orthogonalization, since no MO-coefficients are available
        arma::mat u;
        arma::vec e;
        eig_sym(e, u, m_moldata->s);

        m_moldata->c_fb.alpha() = u * arma::diagmat(1/sqrt(e)) * u.t();
    }

    // Multipole matrices
    {
        // Read the origins used in the operator definitions
        DataSet Set = m_file.openDataSet("MLTPL_ORIG");
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
        arma::mat mp_orig(buf, 3, 3);

        arma::vec dorig = mp_orig.col(1), qorig = mp_orig.col(2);
        if (norm(dorig)!=0.) {
            throw libwfa_exception(k_clazz,
                method, __FILE__, __LINE__, "Dipole moment integrals not centered at origin");
        }

        // Read the operator matrices
        m_moldata->mom.set(0, 0) = m_moldata->s;
        read_mltpl_mat("AO_MLTPL_X",  0, 1);
        read_mltpl_mat("AO_MLTPL_Y",  1, 1);
        read_mltpl_mat("AO_MLTPL_Z",  2, 1);
        read_mltpl_mat("AO_MLTPL_XX", 0, 2);
        read_mltpl_mat("AO_MLTPL_YY", 1, 2);
        read_mltpl_mat("AO_MLTPL_ZZ", 2, 2);

        // Shift the quadrupole operators to the origin of the coordinate system
        if (accu(qorig%qorig)!=0.) {
            std::cout << "Shifting the quadrupole operators to the origin ..." << std::endl;
            for (size_t icart = 0; icart < 3; icart++) {
                arma::mat shift = 2*qorig.at(icart)*m_moldata->mom.get(icart,1) -
                    qorig.at(icart)*qorig.at(icart)*m_moldata->mom.get(0,0);
                m_moldata->mom.set(icart,2) += shift;
            }
        }
    }
}

void molcas_wf_analysis_data::read_input(char *inp) {
    m_input = std::unique_ptr<input_data>(new input_data());

/*  const char* Project = std::getenv("Project");
    std::string inpname(Project);
    inpname.append(".Wfa.Input");
    std::ifstream infile(inpname.c_str());*/

    std::stringstream infile(inp);

    std::cout << "Parsing input file ..." << std::endl;

    bool inwfa = false;
    bool fread;
    std::string str, str4;
    int nread;
    int ctnum_mode = -1;

    while (infile >> str) {
        str4 = str.substr(0,4);
        std::transform(str4.begin(), str4.end(), str4.begin(), ::toupper);

        if (str4=="&WFA") {
            inwfa = true;
        }
        else if (inwfa) {
            if (str4=="H5FI") {
                infile >> str;
                m_input->file_name = str;
            }
            else if (str4=="REFS") {
                infile >> str;
                m_input->refstate = atoi(str.c_str()) - 1;
            }
            else if (str4=="WFAL") {
                infile >> str;
                m_input->wfalevel = atoi(str.c_str());
            }
            else if (str4=="CTNU") {
                infile >> str;
                ctnum_mode = atoi(str.c_str());
            }
            else if (str4=="PROP") {
                m_input->prop_list.clear();

                fread = false;
                while(infile >> str) {
                    if (str=="*") {
                        fread = true;
                        break;
                    }
                    m_input->prop_list.push_back(str);
                }
                if (!fread)
                    throw libwfa_exception(k_clazz, "read_input (PROPLIST)",
                    __FILE__, __LINE__, "Use * to finish input!");
                 m_input->ctnum = true;
            }
            else if (str4=="ATLI") {
                infile >> str;
                nread = atoi(str.c_str());
                m_input->at_lists = std::vector<std::vector<int>>(nread);

                for (int i = 0; i < nread; i++) {
                    while(infile >> str) {
                        if (str=="*") {
                            fread = true;
                            break;
                        }
                        m_input->at_lists[i].push_back(atoi(str.c_str()));
                    }
                    if (!fread)
                        throw libwfa_exception(k_clazz, "read_input (ATLISTS)",
                        __FILE__, __LINE__, "Use * to finish input!");
                     m_input->ctnum = true;
                }
                std::cout << "ATLISTS parsed:" << std::endl;
                for (int i = 0; i < nread; i++) {
                    std::cout << "[ ";
                    for (size_t j = 0; j < m_input->at_lists[i].size(); j++) {
                        std::cout << m_input->at_lists[i][j] << " ";
                    }
                    std::cout << "], ";
                }
                std::cout << std::endl;
            }
            else if (str4=="MULL") m_input->mulliken = true;
            else if (str4=="LOWD") m_input->lowdin = true;
            else if (str4=="NXO")  m_input->nxo = true;
            else if (str4=="EXCI") m_input->exciton = true;
            else if (str4=="DOCT") m_input->ctnum = true;
            else if (str4=="H5OR") m_input->h5orbs = true;

            else if (str4=="DEBU") m_input->debug = true;
            else if (str4=="ADDI") m_input->add_info = true;

            else if (str4=="END") break;

            else {
                std::cout << std::endl << "Unknown keyword: " << str4 << std::endl << std::endl;
                throw libwfa_exception(k_clazz, "read_input", __FILE__, __LINE__, "Unknown keyword");
            }
        }
    }

    // Automatically activate options according to wfalevel
    if (m_input->wfalevel >= 1) {
        if (!m_input->mulliken)
            m_input->lowdin = true;
        m_input->nxo = true;
    }
    if (m_input->wfalevel >= 2) {
        m_input->exciton = true;
    }
    if (m_input->wfalevel >= 3) {
        m_input->ctnum = true;
        m_input->h5orbs = true;
    }
    if (m_input->wfalevel >= 4) {
        m_input->lowdin = true;
        m_input->mulliken = true;
    }

    // Activate fragment-based analysis if fragments are defined and ctnum_mode is given
    if (m_input->at_lists.size() == 0 || ctnum_mode==0) {
        m_input->prop_list = {"Om"};
    }
    else if (ctnum_mode==1) {
        m_input->prop_list = {"Om", "POS", "PR", "DEL", "CT", "CTnt"};
    }
    else if (ctnum_mode==2)  {
        m_input->prop_list = {"Om", "POS", "POSi", "POSf", "PR", "PRi", "PRf",
        "DEL", "COH", "CT", "CTnt"};
    }
    else if (ctnum_mode==3) {
        m_input->prop_list = {"Om", "POSi", "POSf", "PR", "CT",
        "MC", "LC", "MLCT", "LMCT", "LLCT"};
    }
}

std::string molcas_wf_analysis_data::get_mo_types(const H5std_string &setname) {
    DataSet Set = m_file.openDataSet(setname);
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

arma::mat molcas_wf_analysis_data::get_mo_vectors(const H5std_string &setname) {
    DataSet Set = m_file.openDataSet(setname);
    DataSpace Space = Set.getSpace();
    int rank = Space.getSimpleExtentNdims();
    if (rank != 1)
        throw libwfa_exception(k_clazz, "get_mo_vectors", __FILE__, __LINE__, "Inconsistent rank");

    hsize_t dim;
    Space.getSimpleExtentDims(&dim, NULL);
    double* buf = new double[dim];
    Set.read(buf, PredType::NATIVE_DOUBLE);

    arma::mat retmat;
    size_t nsym = m_moldata->nbas.size();
    read_ao_mat(buf, dim, retmat, nsym, 0);
    if (m_moldata->nbas.size() > 1)
        retmat = m_moldata->desym * retmat;

    return retmat;
}

void molcas_wf_analysis_data::cleanup() {

    delete m_input.release();
    delete m_moldata.release();

    for (std::vector<pa_data *>::iterator i = m_pa.begin();
            i != m_pa.end(); i++) {
        delete *i; *i = 0;
    }
    m_pa.clear();
    for (std::vector<cta_data *>::iterator i = m_cta.begin();
            i != m_cta.end(); i++) {
        delete *i; *i = 0;
    }
    m_cta.clear();
}

void molcas_wf_analysis_data::setup_h5core() {
    if (m_h5core.get() != 0) return;

    m_h5core = std::unique_ptr<molcas_export_h5orbs>(new molcas_export_h5orbs(m_file, m_moldata->nbas, m_moldata->desym));
}

void molcas_wf_analysis_data::read_ao_mat(const double *buf, const size_t dim, arma::mat &ao_mat, const size_t nsym, const size_t psym) {
    size_t nbas_t = arma::accu(m_moldata->nbas);
    //size_t nsym = m_moldata->nbas.size();

    if (nsym == 1) { // no symmetry
        if (nbas_t * nbas_t != dim)
            throw libwfa_exception(k_clazz, "read_ao_mat", __FILE__, __LINE__, "Inconsistent AO-matrix (no symmetry)");

        ao_mat = arma::mat(buf, nbas_t, nbas_t);
    } // symmetry
    else {
        hsize_t dim_chk = 0;
        for (size_t isym=0; isym<nsym; isym++)
            dim_chk += m_moldata->nbas(isym) * m_moldata->nbas(isym);

        if (dim_chk != dim)
            throw libwfa_exception(k_clazz, "read_ao_mat", __FILE__, __LINE__, "Inconsistent AO-matrix (symmetry)");

        size_t i=0, j=0, jsym;
        const double *buf_ptr = buf;
        ao_mat = arma::zeros(nbas_t, nbas_t);
        for (size_t isym=0; isym<nsym; isym++) {

            jsym = m_moldata->symmult(isym, psym);

            int ibas = m_moldata->nbas(isym);
            int jbas = m_moldata->nbas(jsym);
            if (ibas*jbas==0) {
                i+=ibas;
                continue;
            }

            j = 0;
            for (size_t jjsym=0; jjsym<jsym; jjsym++) {
                j += m_moldata->nbas(jjsym);
            }

            ao_mat.submat(i, j, i+ibas-1, j+jbas-1) = arma::mat(buf_ptr, ibas, jbas);

            i+=ibas;
            buf_ptr += ibas*jbas;

    }
    if (m_input->debug) {
        ao_mat.print("Output from read_ao_mat");
    }
}
}

void molcas_wf_analysis_data::read_mltpl_mat(const H5std_string &setname, const size_t c, const size_t n) {
    DataSet Set = m_file.openDataSet(setname);
    DataSpace Space = Set.getSpace();
    int rank = Space.getSimpleExtentNdims();
    if (rank != 2)
        throw libwfa_exception(k_clazz, "read_mltpl_mat", __FILE__, __LINE__, "Inconsistent rank");

    hsize_t dims[rank];
    Space.getSimpleExtentDims(dims, NULL);
    hsize_t dimt = dims[0] * dims[1];
    double* buf = new double[dimt];
    Set.read(buf, PredType::NATIVE_DOUBLE);

    read_ao_mat(buf, dimt, m_moldata->mom.set(c, n), 1, 0);
    if (m_input->debug) {
        std::cout << "*** " << c << ", " << n << std::endl;
        m_moldata->mom.get(c,n).print("symm");
    }

    if (m_moldata->nbas.size() > 1) {
        //arma::mat tmp = m_moldata->mom.get(c, n);
        m_moldata->mom.set(c, n) = m_moldata->desym * m_moldata->mom.get(c, n) * m_moldata->desym.t();
        if (m_input->debug) {
            m_moldata->mom.get(c,n).print("desymm");
        }
    }
}

molcas_wf_analysis_data *molcas_setup_wf_analysis_data(char *inp) {
    molcas_wf_analysis_data *h = new molcas_wf_analysis_data(inp);


    // Activate population analyses
    if (h->input()->mulliken) {
        h->init_pop_analysis("mulliken");
    }
    if (h->input()->lowdin) {
        h->init_pop_analysis("lowdin");
    }

    if (h->input()->ctnum) {
        if (h->input()->mulliken) {
            h->init_ctnum_analysis("mulliken");
        }
        if (h->input()->lowdin) {
            h->init_ctnum_analysis("lowdin");
        }
    }

    if (h->input()->nxo) {
        h->activate(molcas_wf_analysis_data::NO);
        h->activate(molcas_wf_analysis_data::NDO);
        h->activate(molcas_wf_analysis_data::NTO);
        h->activate(molcas_wf_analysis_data::SA_NTO);

        size_t norb = 3;
        double thresh = 1.e-5;

        // Setup parameters for orbital print-out; currently use the same parameters
        // for NOs, NDOs, and NTOs
        h->set_orbital_params(orbital_type::NO,  norb, thresh);
        h->set_orbital_params(orbital_type::NDO, norb, thresh);
        h->set_orbital_params(orbital_type::NTO, norb, thresh);

        // Initialize orbital export
        orbital_type::flag_t ot;
        ot.set(orbital_type::no);
        ot.set(orbital_type::ndo);
        ot.set(orbital_type::nto);
        if (h->input()->h5orbs)
            h->init_orbital_export("h5", ot);
    }

    if (h->is_active(wf_analysis_data_i::NTO))
        h->activate(molcas_wf_analysis_data::FORM_EH);
    if (h->is_active(wf_analysis_data_i::NDO))
        h->activate(molcas_wf_analysis_data::FORM_AD);

    if (h->input()->exciton) {
        h->activate(molcas_wf_analysis_data::DENS_MOM);
        h->activate(molcas_wf_analysis_data::EXCITON);
        h->activate(molcas_wf_analysis_data::EXCITON_AD);
    }

    return h;
}

} // namespace libwfa
