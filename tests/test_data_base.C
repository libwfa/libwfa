#include <fstream>
#include "test_data_base.h"

namespace libwfa {

using namespace arma;


bool test_data_base::read_double(const char *testname,
    const char *fname, double &d) {

    std::ifstream in(make_filename(fname).c_str());
    if (in.fail()) return false;

    in >> d;
    if (in.fail()) return false;

    return true;
}


bool test_data_base::read_ab_matrix(const char *testname,
	const char *name, ab_matrix &m) {

	std::string fn(name);
	fn += "_a";
	bool ok = read_matrix(testname, fn.c_str(), m.alpha());

	if (! m.is_alpha_eq_beta()) {
		fn = std::string(name) + std::string("_b");
		ok = ok && read_matrix(testname, fn.c_str(), m.beta());
	}
	return ok;
}


std::string test_data_base::make_filename(const char *fname) const {

    // TODO: make this system independent
    static const char k_sep = '/';

    return std::string(m_prefix + k_sep + fname);
}


} // namespace libwfa
