#include "molcas_export_h5orbs.h"

namespace libwfa {
    
using namespace H5;    
void molcas_export_h5orbs::perform(const std::string &name,
    const arma::mat &c, const arma::vec &e, size_t nocc) {
    
    H5std_string oname(name + "_VECTORS");
    std::transform(oname.begin(), oname.end(), oname.begin(), ::toupper);
    
    std::cout << "\n  *** Exporting " << oname << std::endl;
    
    const int rank = 1;
    hsize_t dim = c.size();
    DataSpace Space(rank, &dim);

    const double *dataptr = c.memptr();

    H5File wfile("molcas.scf.h5", H5F_ACC_RDWR);

    DataSet Set;
    try {
        Set = m_file.createDataSet(oname, PredType::NATIVE_DOUBLE, Space);
    }
    catch( FileIException error ) { // open and overwrite the dataset if it already exists
        std::cout << std::endl << "Overwriting existing dataset " << oname << std::endl;
        Set = m_file.openDataSet(oname);
    }
    Set.write(dataptr, PredType::NATIVE_DOUBLE);
}

void molcas_export_h5orbs::perform(const std::string &name,
    const arma::mat &c_a, const arma::vec &e_a, size_t nocc_a,
    const arma::mat &c_b, const arma::vec &e_b, size_t nocc_b) {
    
    std::cout << "*** Export for unrestricted orbitals not implemented yet!" << std::endl;
}
}