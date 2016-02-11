#include <libwfa/libwfa.h>
#include "molcas_wf_analysis_data.h"

namespace libwfa {
using namespace H5;

const char molcas_wf_analysis_data::k_clazz[] = "molcas_wf_analysis_data";
    
molcas_wf_analysis_data::molcas_wf_analysis_data(H5File file) :
    m_file(file), m_export_dens(EXPORT_NONE), m_export_orbs(EXPORT_NONE) {

    initialize();
    std::cout << "Initialization finished." << std::endl;
}

void molcas_wf_analysis_data::initialize() {
    
//    static const char m_clazz[] = "molcas_wf_analysis_data";    
    static const char method[] = "initialize()";    
    
    Group Grp_main = m_file.openGroup("/");
    hsize_t natoms;
    
    // Read basic info
    {
        Attribute Att_nbas = Grp_main.openAttribute("NBAS");
        DataSpace Spac_nbas = Att_nbas.getSpace();
        hsize_t nnbas;
        Spac_nbas.getSimpleExtentDims(&nnbas, NULL);
        
        int nbas[nnbas];
        Att_nbas.read(PredType::NATIVE_INT, nbas);
    
        int nbas_t = 0;
        std::cout << "nbas:";
        for (int ibas=0; ibas<nnbas; ibas++) {
            std::cout << " " << nbas[ibas];
            nbas_t += nbas[ibas];
        }
        std::cout << std::endl;    
    
        m_moldata = std::auto_ptr<base_data>(new base_data(nbas_t, 0, 0, false));
    }
    
    // Read atom names
    // TODO
    
    // Read atomic numbers/charges. Different in the case of ECPs(?)
    {
        DataSet Set = m_file.openDataSet("CENTER_CHARGES");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 1) {            
            throw libwfa_exception(k_clazz,
                method, __FILE__, __LINE__, "Inconsistent rank");
        }
        
        Space.getSimpleExtentDims(&natoms, NULL);
        std::cout << "-> Number of atoms: " << natoms << std::endl;
        double t[natoms];
        int num_buf[natoms];
        Set.read(&num_buf, PredType::NATIVE_INT);
        
        m_moldata->atoms = std::vector<std::string>(natoms);
        m_moldata->atomic_numbers = arma::uvec(natoms);
        m_moldata->atomic_charges = arma::vec(natoms);
        for (size_t i = 0; i < natoms; i++) {
            m_moldata->atomic_numbers.at(i) = num_buf[i];
            m_moldata->atomic_charges.at(i) = num_buf[i];
            m_moldata->atoms.at(i) = "X";
        }
        std::cout << std::endl;
    }
    
    // Read atomic coordinates
    {
        DataSet Set = m_file.openDataSet("CENTER_COORDINATES");
        DataSpace Space = Set.getSpace();
        int rank = Space.getSimpleExtentNdims();
        if (rank != 2) {            
            throw libwfa_exception(k_clazz,
                method, __FILE__, __LINE__, "Inconsistent rank");
        }
        
        hsize_t dims[rank];
        Space.getSimpleExtentDims(dims, NULL);
        std::cout << "dims: " << dims[0] << " " << dims[1] << std::endl;
        int dim_t = dims[0] * dims[1];
        double coor_buf[dim_t];
        Set.read(&coor_buf, PredType::NATIVE_DOUBLE);
        m_moldata->coordinates = arma::mat(coor_buf, 3, natoms);
        m_moldata->coordinates.print();
    }
    
}

molcas_wf_analysis_data *molcas_setup_wf_analysis_data(H5::H5File file) {
    std::cout << "Starting setup..." << std::endl;
    
    molcas_wf_analysis_data *h = new molcas_wf_analysis_data(file);
    
    return h;
}
} // namespace libwfa