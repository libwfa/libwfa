#include "export_cube_base.h"
#include "libwfa_exception.h"

namespace libwfa {


namespace cube {


grid3d::grid3d() {
    for(size_t i = 0; i < 6; i++) ranges[i] = 0.0;
    for(size_t i = 0; i < 3; i++) npts[i] = 0;
}


void grid3d::set_xrange(double min, double max, unsigned int npts_) {

    if (min >= max) {
        throw libwfa_exception("cube::grid3d",
                "set_xrange(double, double, unsigned int)",
                __FILE__, __LINE__, "min > max");
    }
    if (npts_ == 0) {
        throw libwfa_exception("cube::grid3d",
                "set_xrange(double, double, unsigned int)",
                __FILE__, __LINE__, "npts");
    }

    ranges[0] = min;
    ranges[1] = max;
    npts[0] = npts_;
}


void grid3d::set_yrange(double min, double max, unsigned int npts_) {

    if (min >= max) {
        throw libwfa_exception("cube::grid3d",
                "set_xrange(double, double, unsigned int)",
                __FILE__, __LINE__, "min > max");
    }
    if (npts_ == 0) {
        throw libwfa_exception("cube::grid3d",
                "set_xrange(double, double, unsigned int)",
                __FILE__, __LINE__, "npts");
    }

    ranges[2] = min;
    ranges[3] = max;
    npts[1] = npts_;
}


void grid3d::set_zrange(double min, double max, unsigned int npts_) {

    if (min >= max) {
        throw libwfa_exception("cube::grid3d",
                "set_xrange(double, double, unsigned int)",
                __FILE__, __LINE__, "min > max");
    }
    if (npts_ == 0) {
        throw libwfa_exception("cube::grid3d",
                "set_xrange(double, double, unsigned int)",
                __FILE__, __LINE__, "npts");
    }

    ranges[4] = min;
    ranges[5] = max;
    npts[2] = npts_;
}


void grid3d::check() const {

    char dim = 'x';
    for (size_t i = 0; i < 3; i++, dim++) {
        if (ranges[2 * i] >= ranges[2 * i + 1]) {
            std::ostringstream ss; ss << dim << ": min < max";
            throw libwfa_exception("grid3d", "check()",
                     __FILE__, __LINE__, ss.str().c_str());
        }
        if (npts[i] == 0) {
            std::ostringstream ss; ss << dim << ": npts";
            throw libwfa_exception("grid3d", "check()",
                     __FILE__, __LINE__, ss.str().c_str());
        }
    }
}


const data_type data_type::orb_a = data_type(+1);
const data_type data_type::orb_b = data_type(-1);
const data_type data_type::sdm_a = data_type(+2);
const data_type data_type::sdm_b = data_type(-2);
const data_type data_type::tdm_a = data_type(+3);
const data_type data_type::tdm_b = data_type(-3);
const data_type data_type::adm_a = data_type(+4);
const data_type data_type::adm_b = data_type(-4);
const data_type data_type::ddm_a = data_type(+5);
const data_type data_type::ddm_b = data_type(-5);


} // namespace cube


} // namespace libwfa



