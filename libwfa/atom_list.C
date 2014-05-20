#include <iomanip>
#include "atom_list.h"
#include "libwfa_exception.h"


namespace libwfa {


atom_list::atom_list(size_t natoms) :
    m_natoms(natoms), m_atnum(0), m_coords(0), m_alloc(true) {

    m_atnum = new int[m_natoms];
    m_coords = new double[m_natoms * 3];
}


atom_list::atom_list(size_t natoms, const int *atnum, const double *coords) :
    m_natoms(natoms), m_atnum(0), m_coords(0), m_alloc(true) {

    m_atnum = new int[m_natoms];
    m_coords = new double[m_natoms * 3];

    memcpy(m_atnum, atnum, m_natoms * sizeof(int));
    memcpy(m_coords, coords, m_natoms * sizeof(double) * 3);
}


atom_list::atom_list(size_t natoms, int *atnum, double *coords, bool do_copy) :
    m_natoms(natoms), m_atnum(atnum), m_coords(coords), m_alloc(do_copy) {

}


void atom_list::set_atom(size_t idx, int atnum, const double (&x)[3]) {

#ifdef LIBWFA_DEBUG
    check_idx(idx, m_natoms);
#endif

    m_atnum[idx] = atnum;
    for (size_t i = 0, j = idx * 3; i < 3; i++, j++) m_coords[j] = x[i];
}


void atom_list::check_idx(size_t idx, size_t limit) {

    if (idx >= limit) {
        throw libwfa_exception("atom_list", "check_idx()",
                    __FILE__, __LINE__, "idx");
     }
}


std::ostream &operator<<(std::ostream &out, const atom_list &lst) {

    for (size_t i = 0; i < lst.natoms(); i++) {
        out << std::setw(5) << lst.atomic_number(i);
        out << std::setw(12) << std::setprecision(6) << std::fixed;
        const double *pos = lst.position(i);
        for (size_t j = 0; j < 3; j++) out << " " << pos[j];
        out << std::endl;
    }

    return out;
}


} // namespace libwfa

