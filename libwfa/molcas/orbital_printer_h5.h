#ifndef LIBWFA_ORBITAL_PRINTER_H5_H
#define LIBWFA_ORBITAL_PRINTER_H5_H

#include "H5Cpp.h"
#include <libwfa/export/orbital_printer_i.h>

namespace libwfa {

class orbital_printer_h5 : public orbital_printer_i {
public:
    static const char k_clazz[]; //!< Class name
    
public:
    typedef orbital_type::flag_t ot_flag; //!< Flag for orbital types};
    
private:
    H5::H5File m_file; //!< HDF5 file
    std::string m_id; //!< ID / name of orbitals
    ot_flag m_ot; //!< Flag which orbital types are printed

public:
    /** \brief Constructor
        \param xx HDF5 file
        \param id ID / name of orbitals
        \param ot Flag which orbital types to export
     **/
    orbital_printer_h5(H5::H5File &file, const std::string &id, const ot_flag &ot = ot_flag(orbital_type::ALL)) :
        m_file(file), m_id(id), m_ot(ot) {}
    
    /** \brief Destructor
     **/
    virtual ~orbital_printer_h5() { }    

    /** \copydoc orbital_printer_i::perform(orbital_type, const orbital_data &,
            const orbital_selector &)
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb, const orbital_selector &s);

    /** \copydoc orbital_printer_i::perform(orbital_type, const orbital_data &,
            const orbital_selector &, const orbital_data &,
            const orbital_selector &)
     **/
    virtual void perform(orbital_type type,
            const orbital_data &orb_a, const orbital_selector &s_a,
            const orbital_data &orb_b, const orbital_selector &s_b);
};

} // namespace libwfa

#endif // LIBWFA_ORBITAL_PRINTER_H5_H    