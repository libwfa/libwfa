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

    size_t npts, batchsz = 100, nx[3];
    double *gpts = new double[batchsz * 3];
    for (size_t ipt = 0; ipt != npts; ipt += batchsz) {

        size_t sz = std::min(batchsz, npts - ipt + 1);
        size_t j = 1, pos = 3, d = 2;
        while (j != sz) {
            for (size_t k = 0, p0 = pos - 3, p1 = pos; k != 3; k++, p0++, p1++) {
                gpts[p1] = gpts[p0] + m_grid.direction(d, k);
            }
            j++; i++
            if (j > nx[1]) {
                d += nx[i]
            }
            else if (j > nx[2]) {
            }
        }
        for (size_t j = 0; j != sz; j++) {
            size_t d = 2;

        }
    }

    // Loop over batches

    // Generate grid points for batch

    // Evaluate basis functions on grid points

    // batch size * nbasis
    //
    EvlBaD();

    // Compute density on grid

    // batch == output
    // jChi
    // p == density data
    // n == basis
    // npper == batch size
    // nsets == number of densities

    DnMesh();
}



void export_cube_qchem::perform(const std::string &name,
        const std::vector<size_t> &idx, const arma::Mat<double> &vecs) {

    throw libwfa_exception("export_cube_qchem", "perform()",
             __FILE__, __LINE__, "NIY");

}


} // namespace libwfa



