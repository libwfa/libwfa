#include "multipol_con.h"

namespace libwfa{

using namespace arma;


double multipol_con::perform(const Mat<double> &dm,
		const std::string &op1, const std::string &op2) const{

	const Mat<double> &m1 = retrieve_op(op1);
	const Mat<double> &m2 = retrieve_op(op2);

	return accu((dm * m1) % (m2 * dm));

}//end fct


double multipol_con::perform (const Mat<double> &dm, const std::string &op) const{

	const Mat<double> &m = retrieve_op(op);

	return trace(dm * m);
} //end fct


const Mat<double> &multipol_con::retrieve_op(const std::string &op) const{

	if (op == "x") return m_x;
	else if (op == "y") return m_y;
	else if (op == "z") return m_z;
	else if (op == "xx") return m_xx;
	else if (op == "yy") return m_yy;
	else if (op == "zz") return m_zz;
	else if (op == "s") return m_s;
	else
		throw 1;
}//end fct









}//end namespace libwfa
