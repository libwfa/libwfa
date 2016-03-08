#include "molcas_export_h5orbs.h"

namespace libwfa {

using namespace H5;
void molcas_export_h5orbs::perform(const std::string &name,
    const arma::mat &c, const arma::vec &e, size_t nocc) {

    // Export coefficients
    h5_vec("DESYM_" + name + "_VECTORS", vectorise(c));

    // Export eigenvalues - these are interpreted as occupations
    h5_vec("DESYM_" + name + "_OCCUPATIONS", e);
}

void molcas_export_h5orbs::perform(const std::string &name,
    const arma::mat &c_a, const arma::vec &e_a, size_t nocc_a,
    const arma::mat &c_b, const arma::vec &e_b, size_t nocc_b) {

    // Export alpha coefficients
    h5_vec("DESYM_" + name + "_ALPHA_VECTORS", vectorise(c_a));

    // Export eigenvalues - these are interpreted as occupations
    h5_vec("DESYM_" + name + "_ALPHA_OCCUPATIONS", e_a);

    // Export beta coefficients
    h5_vec("DESYM_" + name + "_BETA_VECTORS", vectorise(c_b));

    // Export eigenvalues - these are interpreted as occupations
    h5_vec("DESYM_" + name + "_BETA_OCCUPATIONS", e_b);
}

void molcas_export_h5orbs::h5_vec(H5std_string oname, const arma::vec data) {
    std::transform(oname.begin(), oname.end(), oname.begin(), ::toupper);

    const hsize_t dim = data.size();
    const double *dataptr = data.memptr();
    const hsize_t rank = 1;
    DataSet Set;

    try {
        DataSpace Space(rank, &dim);
        Set = m_file.createDataSet(oname, PredType::NATIVE_DOUBLE, Space);
    }
    catch( FileIException error ) { // open and overwrite the dataset if it already exists
        std::cout << std::endl << "Overwriting existing dataset " << oname << std::endl;
        Set = m_file.openDataSet(oname);
    }
    Set.write(dataptr, PredType::NATIVE_DOUBLE);
}

} // namespace libwfa