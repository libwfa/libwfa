#include "export_cube_base.h"
#include "libwfa_exception.h"

namespace libwfa {


void export_cube_qchem::perform(const std::string &name,
        const arma::Mat<double> &mat) {

    throw libwfa_exception("export_cube_qchem", "perform()",
             __FILE__, __LINE__, "NIY");

    // - get_carts to retrieve the molecular coordinates
    // - pntper to determine the batching size
    // - MkAtmO shell offsets per atom
    // - SgS1Gd significant shells on the current grid
}


void export_cube_qchem::perform(const std::string &name,
        const std::vector<size_t> &idx, const arma::Mat<double> &vecs) {

    throw libwfa_exception("export_cube_qchem", "perform()",
             __FILE__, __LINE__, "NIY");

}


} // namespace libwfa



