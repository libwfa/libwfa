#include <sstream>
#include "state_info.h"

namespace libwfa {


const spin spin::unspecified = spin(0);
const spin spin::singlet = spin(1);
const spin spin::doublet = spin(2);
const spin spin::triplet = spin(3);
const spin spin::quartet = spin(4);
const spin spin::quintet = spin(5);


std::string state_info::convert(char sep) const {

    std::ostringstream ss;
    ss << stateno;
    if (multiplicity != spin::unspecified) {
        ss << sep << multiplicity.multiplicity();
    }
    ss << sep << irrep;

    return ss.str();
}


} // namespace libwfa


