#include "dm_type.h"

namespace libwfa {


const dm_type dm_type::sdm = dm_type(0);
const dm_type dm_type::tdm = dm_type(1);
const dm_type dm_type::adm = dm_type(2);
const dm_type dm_type::ddm = dm_type(3);
const dm_type dm_type::edm = dm_type(4);
const dm_type dm_type::hdm = dm_type(5);


std::string dm_type::convert() const {
    switch (m_type) {
    case 0: return "sdm";
    case 1: return "tdm";
    case 2: return "adm";
    case 3: return "ddm";
    case 4: return "edm";
    case 5: return "hdm";
    default: return "unknown";
    }
}


} // namespace libwfa

