#include <libwfa/libwfa_exception.h>
#include "orbital_printer_h5.h"

namespace libwfa {

using namespace H5;

const char orbital_printer_h5::k_clazz[] = "orbital_printer_h5";

void orbital_printer_h5::perform(orbital_type type,
        const orbital_data &orb, const orbital_selector &s) {
    std::cout << "\n  *** Exporting restricted orbitals: " << m_id << std::endl;
    
    // TOOD:desymmetrize
    
    const int rank = 1;
    hsize_t dim = 5;
    DataSpace Space(rank, &dim);
    
    double vals[] = {1., 2., 3., 4., 5.};
    
    std::cout << "fptmp A" << std::endl;
    
    H5File wfile("molcas.scf.h5", H5F_ACC_RDWR);
    
    // TODO: Check if dataset already exists ...
    DataSet Set(m_file.createDataSet("a2", PredType::NATIVE_DOUBLE, Space));
    
    Set.write(vals, PredType::NATIVE_DOUBLE);

    std::cout << "fptmp B" << std::endl;
}

void orbital_printer_h5::perform(orbital_type type,
        const orbital_data &orb_a, const orbital_selector &s_a,
        const orbital_data &orb_b, const orbital_selector &s_b) {
    std::cout << " *** Exporting unrestricted orbitals" << std::endl;
}
} // namespace libwfa