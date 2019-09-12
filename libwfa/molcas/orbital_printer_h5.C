#include <libwfa/libwfa_exception.h>
#include "orbital_printer_h5.h"

namespace libwfa {

using namespace H5;
using namespace arma;

const char orbital_printer_h5::k_clazz[] = "orbital_printer_h5";

void orbital_printer_h5::perform(orbital_type type,
        const orbital_data &orb, const orbital_selector &s) {

    // TOOD:desymmetrize

    H5std_string name(m_id + "_" + type.convert_upper() + "_VECTORS");
    std::cout << "\n  *** Exporting " << type << ": " << name << std::endl;

    //size_t no = s.n_selected(true);

    //uvec idx = s.get_selected_arma();
    //if (idx.n_rows != 0) {
        mat c = orb.get_coeff();//.cols(idx);

        const int rank = 1;
        hsize_t dim = c.size();
        DataSpace Space(rank, &dim);

        double *dataptr = c.memptr();

        std::cout << "fptmp A" << std::endl;

        H5File wfile("molcas.scf.h5", H5F_ACC_RDWR);

        DataSet Set;
        try {
            Set = m_file.createDataSet(name, PredType::NATIVE_DOUBLE, Space);
        }
        catch( FileIException error ) { // open and overwrite the dataset if it already exists
            std::cout << std::endl << "Overwriting existing dataset " << name << std::endl;
            Set = m_file.openDataSet(name);
        }
        Set.write(dataptr, PredType::NATIVE_DOUBLE);

        std::cout << "fptmp B" << std::endl;
    //}
}

void orbital_printer_h5::perform(orbital_type type,
        const orbital_data &orb_a, const orbital_selector &s_a,
        const orbital_data &orb_b, const orbital_selector &s_b) {
    std::cout << " *** Exporting unrestricted orbitals" << std::endl;
}
} // namespace libwfa