#ifndef EX_ANALYSE_H_
#define EX_ANALYSE_H_

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/multipol_con_i.h>


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
    create an object for the calculation and the class,
    call the perform function once, then get your desired results
    from the corresponding arrays.
    Access to the arrays is provided by get_ functions..

    \ingroup libwfa
 **/
class ex_analyse {
private:
    double rh[3][2]; //!< \f$<r_h>\f$ for x,y,z and alpha/beta
    double rh2[3][2];//!< \f$<(r_h)²>\f$ for x,y,z and alpha/beta
    double re[3][2];//!< \f$<r_e>\f$ for x,y,z and alpha/beta
    double re2[3][3];//!< \f$<(r_e)²>\f$ for x,y,z and alpha/beta
    double rhre[3][2];//!< \f$<r_e*r_h>\f$ for x,y,z and alpha/beta
	double sep[2];//!< \f$|<r_h - r_e>|\f$ for alpha/beta
	double dex_c[3][2];//!< \f$<(r_h)²>-2<r_h*r_e>+<(r_e)²>\f$ for x,y,z and a,b
	double dex_tot[2];//!< \f$sqrt(<d_ex,x>²+<d_ex,y>²+<d_ex,z>)\f$ for a, b
	double sig_h[2];//!< \f$sqrt(<(r_h)²>-<r_h>²)\f$ for alpha/beta
	double sig_e[2];//!< \f$sqrt(<(r_e)²>-<r_e>²)\f$ for alpha/beta
	double cov[2];//!< \f$<r_h*r_e>-<r_h>*<r_e>\f$ for alpha/beta
	double corr[2];//!< \f$cov/((\omega_h)*(\omega)_e)\f$ for alpha/beta
	bool m_aeqb;//!< is alpha equal to beta

public:
	/** \brief Forms all needed expected values over all spacial coordinates
         and spins.
    \param tdm Transition Density Matrix
    \param om Omega Omega Matrix
    \param name Reference to the Contract class interface
	 **/
	void ex_form(const ab_matrix &tdm, const ab_matrix &om,
			const multipol_con_i &name);
	/** \brief Calculates distance \f$<(r_h)²>-2<r_h*r_e>+<(r_e)²>\f$ over one
         spacial coordinate and for a specific spin.
         \param koord Char to choose the sp. coord.--> x,y,z
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_d_ex_c(char koord, char spin);
	/** \brief Calculates the total distance
	      \f$sqrt(<d_ex,x>²+<d_ex,y>²+<d_ex,z>)\f$ for a specific spin.
         \param koord Char to choose the sp. coord.--> x,y,z
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_d_ex_tot(char spin);
	/** \brief Calculates the separation \f$|<r_h - r_e>|\f$.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_mean_sep(char spin);
	/** \brief Calculates \f$\omega=sqrt(<(r_h)²>-<r_h>²)\f$ for hole.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_sig_h(char spin);
	/** \brief Calculates \f$\omega=sqrt(<(r_e)²>-<r_e>²)\f$ for electrons.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_sig_e(char spin);
	/** \brief Calculates the covariance\f$<r_h*r_e>-<r_h>*<r_e>\f$.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_cov(char spin);
	/** \brief Calculates the correlation \f$cov/((\omega_h)*(\omega)_e)\f$.
         \param spin Char to choose the spin --> a,b
         \return result
	 **/
	double ex_corr(char spin);

public:
	/** \brief Returns the exp. value of <r_h> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_rh(char koord, char spin);
	/** \brief Returns the exp. value of <r_e> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_re(char koord, char spin);
	/** \brief Returns the exp. value of <(r_h)²> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_rh2(char koord, char spin);
	/** \brief Returns the exp. value of <(r_e)²> for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_re2(char koord, char spin);
	/** \brief Returns the exp. value of <r_h*r_e> for a sp. coordinate and spin.
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
			const multipol_con_i &name);
	/** \brief Returns the separation for a sp. spin.
        \param spin Char to choose the spin --> a,b
	 **/
	double get_sep(char spin);
	/** \brief Returns the distance for a sp. coordinate and spin.
        \param koord Char to choose the sp. coord.--> x,y,z
        \param spin Char to choose the spin --> a,b
	 **/
	double get_dex_c (char coord, char spin);
	/** \brief Returns the over all distance for a sp. spin.
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
