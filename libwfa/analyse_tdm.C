#include "analyse_tdm.h"
#include "transformations_dm.h"

namespace libwfa {

using namespace arma;


void analyse_tdm::perform(const ab_matrix &tdm, ab_matrix_pair &av,
        export_densities_i &dm_print, export_orbitals_i &nto_print,
        nto_data_i &prn, ctnum_data_i &prct) const {

    dm_print.perform(dm_type::tdm, tdm);

    m_nto.perform(tdm, av, dm_print, nto_print, prn);

    m_ct.perform(tdm, prct);
}


} // end namespace
