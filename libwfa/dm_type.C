#include "dm_type.h"

namespace libwfa {


const dm_type dm_type::state = dm_type(0);
const dm_type dm_type::transition = dm_type(1);
const dm_type dm_type::difference = dm_type(2);
const dm_type dm_type::attach = dm_type(3);
const dm_type dm_type::detach = dm_type(4);
const dm_type dm_type::particle = dm_type(5);
const dm_type dm_type::hole = dm_type(6);


std::string dm_type::convert() const {
    switch (m_type) {
    case 0: return "sdm";
    case 1: return "tdm";
    case 2: return "ddm";
    case 3: return "att";
    case 4: return "det";
    case 5: return "edm";
    case 6: return "hdm";
    default: return "unknown";
    }
}


} // namespace libwfa

