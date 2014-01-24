#include "density_type.h"

namespace libwfa {


const density_type density_type::state = density_type(0);
const density_type density_type::transition = density_type(1);
const density_type density_type::difference = density_type(2);
const density_type density_type::attach = density_type(3);
const density_type density_type::detach = density_type(4);
const density_type density_type::particle = density_type(5);
const density_type density_type::hole = density_type(6);


std::string density_type::convert() const {
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

