#include "orbital_type.h"

namespace libwfa {


const orbital_type orbital_type::mo = orbital_type(0);
const orbital_type orbital_type::no = orbital_type(1);
const orbital_type orbital_type::nto = orbital_type(2);
const orbital_type orbital_type::ndo = orbital_type(3);


std::string orbital_type::convert() const {
    switch (m_type) {
    case 0: return "mo";
    case 1: return "no";
    case 2: return "nto";
    case 3: return "ndo";
    default: return "unknown";
    }
}


std::ostream &operator<<(std::ostream &out, orbital_type type) {

    if (type == orbital_type::mo)
        out << "Molecular orbitals";
    else if (type == orbital_type::no)
        out << "Natural orbitals";
    else if (type == orbital_type::nto)
        out << "Natural transition orbitals";
    else if (type == orbital_type::ndo)
        out << "Natural difference orbitals";
    else
        out << "Unknown orbitals";
    return out;
}


} // namespace libwfa

