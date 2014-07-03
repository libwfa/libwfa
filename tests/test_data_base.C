#include <fstream>
#include "test_data_base.h"

namespace libwfa {

using namespace arma;


bool test_data_base::read_double(const char *testname,
    const char *fname, double &d) const {

    std::ifstream in(make_filename(fname).c_str());
    if (in.fail()) return false;

    in >> d;
    if (in.fail()) return false;

    return true;
}


std::string test_data_base::make_filename(const char *fname) const {

    // TODO: make this system independent
    static const char k_sep = '/';

    return std::string(m_prefix + k_sep + fname);
}


} // namespace libwfa
