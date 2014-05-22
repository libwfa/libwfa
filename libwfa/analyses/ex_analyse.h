#ifndef EX_ANALYSE_H_
#define EX_ANALYSE_H_

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/ab_vector.h>

namespace libwfa {

/**	\brief Class for methods to analyse excitons TDM.

	This class contains several methods to analyse exciton TDMs.
	You can:
	- calculate the expected value for the multipol operators.

 	\ingroup libwfa
 **/
class ex_analyse {
	private:

	public:
		/** \brief Calculates the expected value for the multipol of an exciton and returns either the alpha or beta value.
			\param[in] tdm Transition Density Matrix
			\param[in] op1 Multipoloperator for hole
			\param[in] op2 Multipoloperator for electron
			\param[in] om Omega Matrix
			\param[in] spin Char Switch to set either 'a'=alpha or 'b'=beta
		 **/
		double ex_multip (const ab_matrix &tdm, const ab_matrix &op1, const ab_matrix &op2, ab_matrix &om, char spin);


}; // end Class

} // end namespace libwfa


#endif
