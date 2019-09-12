#include "molcas_export_h5orbs.h"

namespace libwfa {

using namespace H5;

void molcas_export_h5orbs::perform(const std::string &name,
    const arma::mat &c, const arma::vec &e, size_t nocc) {
    // Export coefficients
    h5_coeff(name, c);

    // Export eigenvalues - these are interpreted as occupations
    h5_occs(name, e);
}

void molcas_export_h5orbs::perform(const std::string &name,
    const arma::mat &c_a, const arma::vec &e_a, size_t nocc_a,
    const arma::mat &c_b, const arma::vec &e_b, size_t nocc_b) {

    if (c_a.size() > 0) {
        // Export alpha coefficients
        h5_coeff(name + "_ALPHA", c_a);

        // Export eigenvalues - these are interpreted as occupations
        h5_occs(name + "_ALPHA", e_a);
    }

    if (c_b.size() > 0) {
        // Export beta coefficients
        h5_coeff(name + "_BETA", c_b);

        // Export eigenvalues - these are interpreted as occupations
        h5_occs(name + "_BETA", e_b);
    }
}

void molcas_export_h5orbs::setup(H5File &file) {
        try {
            m_group = file.createGroup("/WFA");

            H5std_string descr = "Wavefunction analysis data produced by wfa_molcas.x";
            DataSpace Space(H5S_SCALAR);
            StrType strtype(PredType::C_S1, descr.size());

            Attribute Att = m_group.createAttribute("DESCRIPTION", strtype, Space);
            Att.write(strtype, descr);
            std::cout << std::endl << "HDF5 data group WFA created." << std::endl;
        }
        catch( H5::FileIException error ) {
            std::cout << std::endl << "WARNING: Overwriting data in existing group WFA " << std::endl;
            m_group = file.openGroup("/WFA");
        }
}

void molcas_export_h5orbs::h5_vec(H5std_string oname, const H5std_string descr, const arma::vec data) {
    std::transform(oname.begin(), oname.end(), oname.begin(), ::toupper);

    const hsize_t dim = data.size();
    const double *dataptr = data.memptr();
    const hsize_t rank = 1;
    DataSet Set;

    try {

        {
            DataSpace Space(rank, &dim);
            Set = m_group.createDataSet(oname, PredType::NATIVE_DOUBLE, Space);
        }
        { // Create a description as attribute
            DataSpace Space(H5S_SCALAR);
            StrType strtype(PredType::C_S1, descr.size());

            Attribute Att = Set.createAttribute("DESCRIPTION", strtype, Space);
            Att.write(strtype, descr);
        }

    }
    catch( GroupIException error ) { // open and overwrite the dataset if it already exists
        //std::cout << "WARNING: Overwriting existing dataset " << oname << std::endl;
        Set = m_group.openDataSet(oname);
    }

    Set.write(dataptr, PredType::NATIVE_DOUBLE);
}

void molcas_export_h5orbs::h5_coeff(H5std_string oname, const arma::mat c) {
    H5std_string descr = "Desymmetrized orbital coefficients";
    h5_vec("DESYM_" + oname + "_VECTORS", descr, vectorise(c));
}

void molcas_export_h5orbs::h5_occs(H5std_string oname, const arma::vec e) {
    H5std_string descr = "Orbital eigenvalues";
    h5_vec("DESYM_" + oname + "_OCCUPATIONS", descr, e);
}


} // namespace libwfa