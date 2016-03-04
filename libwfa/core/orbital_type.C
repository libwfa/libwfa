#include "orbital_type.h"

namespace libwfa {


const orbital_type orbital_type::mo = orbital_type(0);
const orbital_type orbital_type::no = orbital_type(1);
const orbital_type orbital_type::ndo = orbital_type(2);
const orbital_type orbital_type::nto = orbital_type(3);
const orbital_type orbital_type::dyson = orbital_type(4);


std::string orbital_type::convert() const {
    switch (m_type) {
    case 0: return "mo";
    case 1: return "no";
    case 2: return "ndo";
    case 3: return "nto";
    case 4: return "dyson";
    default: return "unknown";
    }
}

std::string orbital_type::convert_upper() const {
    switch (m_type) {
    case 0: return "MO";
    case 1: return "NO";
    case 2: return "NDO";
    case 3: return "NTO";
    case 4: return "DYSON";
    default: return "UNKNOWN";
    }
}

std::ostream &operator<<(std::ostream &out, orbital_type type) {

    if (type == orbital_type::mo)
        out << "Molecular orbitals";
    else if (type == orbital_type::no)
        out << "Natural orbitals";
    else if (type == orbital_type::ndo)
        out << "Natural difference orbitals";
    else if (type == orbital_type::nto)
        out << "Natural transition orbitals";
    else if (type == orbital_type::ndo)
        out << "Dyson orbitals";
    else
        out << "Unknown orbitals";
    return out;
}


} // namespace libwfa

