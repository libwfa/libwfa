#include <libwfa/libwfa.h>
#include "molcas_wf_analysis_data.h"

namespace libwfa {
using namespace H5;
    
molcas_wf_analysis_data::molcas_wf_analysis_data() :
    m_export_dens(EXPORT_NONE), m_export_orbs(EXPORT_NONE) {

    initialize();
}

void molcas_wf_analysis_data::initialize() {
}

molcas_wf_analysis_data *molcas_setup_wf_analysis_data(H5File file) {
    std::cout << "Starting setup..." << std::endl;
    Group Grp_main = file.openGroup("/");
    
    Attribute Att_nbas = Grp_main.openAttribute("NBAS");
    DataSpace Spac_nbas = Att_nbas.getSpace();
    hsize_t nnbas;
    Spac_nbas.getSimpleExtentDims(&nnbas, NULL);
    std::cout << "nnbas: " << nnbas;
    std::cout << std::endl;
    
    int nbas[nnbas];
    Att_nbas.read(PredType::NATIVE_INT, nbas);
    std::cout << "nbas:";
    for (int ibas=0; ibas<nnbas; ibas++)
        std::cout << " " << nbas[ibas];
    std::cout << std::endl;    
}
} // namespace libwfa