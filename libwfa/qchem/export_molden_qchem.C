#include <libgen/libgen.h>
#include "export_molden_qchem.h"

namespace libwfa {


void export_molden_qchem::perform(const std::string &name,
    const ab_matrix &coeff, const ab_vector &ene,
    size_t nocc_a, size_t nocc_b) {

    std::string filename(name + ".mo");
    FILE *fp = QOpen(filename, "w");
    WriteMoldenATOMS(fp);
    WriteMoldenGTO(fp);

    if (coeff.is_alpha_eq_beta()) {

        size_t nv_a = ene.alpha().n_cols - nocc_a;
        double *c_a = const_cast<double *>(coeff.alpha().mem);
        double *e_a = const_cast<double *>(ene.alpha().mem);
        WriteMoldenMOgen(c_a, 0, e_a, 0,
                nocc_a, nv_a, 0, 0, nocc_a, nv_a, 0, 0);
    }
    else {

        size_t nv_a = ene.alpha().n_cols - nocc_a;
        size_t nv_b = ene.beta().n_cols - nocc_b;
        double *c_a = const_cast<double *>(coeff.alpha().mem);
        double *c_b = const_cast<double *>(coeff.beta().mem);
        double *e_a = const_cast<double *>(ene.alpha().mem);
        double *e_b = const_cast<double *>(ene.beta().mem);
        WriteMoldenMOgen(c_a, c_b, e_a, e_b,
                nocc_a, nv_a, nocc_b, nv_b, nocc_a, nv_a, nocc_b, nv_b);
    }

}

} // namespace libwfa



