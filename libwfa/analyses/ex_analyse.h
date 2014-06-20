#ifndef EX_ANALYSE_H_
#define EX_ANALYSE_H_

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/contract_i.h>


namespace libwfa {
using namespace arma;

/** \brief Class for methods to analyse excitons TDM.

    This class contains several methods to analyse exciton TDMs.
    You can:
    - calls the contract class to calc. the exp. val. and saves them
    - calculate the averaged expected disctance
    - calculate the mean seperation
    - calculate the sigma for either the hole or electron
    - calculate the correlation factor

    To use this class you will have to form the needed matrices,
    create an object for the class, call the perform function once,
    then get your desired results from the corresponding arrays.
    Access to the arrays is provided by get_ functions.

    \ingroup libwfa
 **/
class ex_analyse {
private:
	double rh[3][2], re[3][2], rh2[3][2], re2[3][2], rhre[3][2];
	double sep[2];
	double dex_c[3][2];
	double dex_tot[2];
	double sig_h[2];
	double sig_e[2];
	double cov[2];
	double corr[2];
	bool m_aeqb;

public:
	/** \brief Forms all needed expected values over all spacial coordinates
         and spins.
    \param tdm Transition Density Matrix
    \param om Omega Omega Matrix
    \param name Reference to the Contract class interface
	 **/
	void ex_form(const ab_matrix &tdm, const ab_matrix &om,
			const contract_i &name);
	/** \brief Calculates the expected distance over one
         spacial coordinate and for a specific spin.
         \param koord Char to choose the sp. coord.--> x,y,z
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_d_ex_c(char koord, char spin);
	/** \brief Calculates the expected distance over all
         spacial coordinate and for a specific spin.
         \param koord Char to choose the sp. coord.--> x,y,z
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_d_ex_tot(char spin);
	/** \brief Calculates the <rh>-<re> vector and returns its quantity.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_mean_sep(char spin);
	/** \brief Calculates the sigma for the hole.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_sig_h(char spin);
	/** \brief Calculates the sigma for the electron.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_sig_e(char spin);
	/** \brief Calculates the covalenz of rh, re.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_cov(char spin);
	/** \brief Calculates the correlation of rh, re.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_corr(char spin);

public:
	/** \brief Returns the exp. value of <rh> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_rh(char koord, char spin);
	/** \brief Returns the exp. value of <re> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_re(char koord, char spin);
	/** \brief Returns the exp. value of <rh2> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_rh2(char koord, char spin);
	/** \brief Returns the exp. value of <re2> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_re2(char koord, char spin);
	/** \brief Returns the exp. value of <rhre> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_rhre(char koord, char spin);

	/** \brief Performs every analysis in this class and saves it in arrays.
        \param tdm Transition Density Matrix
        \param om Omega Omega Matrix
        \param name Reference to the Contract class interface
	 **/

	void perform(const ab_matrix &tdm, const ab_matrix &om,
			const contract_i &name);
	/** \brief Returns the separation for a sp. spin.
        \param spin Char to choose the spin --> a,b
	 **/
	double get_sep(char spin);
	/** \brief Returns the distance for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_dex_c (char coord, char spin);
	/** \brief Returns the distance for a sp. spin.
        \param spin Char to choose the spin --> a,b
	 **/
	double get_dex_tot (char spin);
	/** \brief Returns the sigma for the hole for a sp. spin.
            \param spin Char to choose the spin --> a,b
	 **/
	double get_sig_h (char spin);
	/** \brief Returns the sigma for the electron for a sp. spin.
            \param spin Char to choose the spin --> a,b
	 **/
	double get_sig_e (char spin);
	/** \brief Returns the covariance for a sp. spin.
            \param spin Char to choose the spin --> a,b
	 **/
	double get_cov (char spin);
	/** \brief Returns the correlation for a sp. spin.
            \param spin Char to choose the spin --> a,b
	 **/
	double get_corr (char spin);
	/** \brief Returns bool if \f$\alpha == \beta\f$
	 **/
	bool aeqb() { return m_aeqb; }

private:
	/**\brief Returns the coordinates as indices
	   \param coord Char to choose coordinate x,y,z
	 **/
	static size_t determine_coord(char coord);
	/**\brief Returns the spin as indices
	   \param spin Char to choose spin a,b
	 **/
	static size_t determine_spin(char spin);
}; // end Class

} // end namespace libwfa


#endif
