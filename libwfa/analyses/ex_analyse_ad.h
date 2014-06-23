#ifndef EX_ANALYSE_AD_H_
#define EX_ANALYSE_AD_H_

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/multipol_con_i.h>

namespace libwfa {

/** \brief Class for methods to analyse excitons through attach/detach matrices.

    This class contains several methods to analyse exciton TDMs.
    You can:
    - calculate the expected value for the multipol operators and saves them
    - calculate the mean seperation
    - calculate the sigma for either the hole or electron

    To use this class you will have to form the needed matrices,
    create an object for the calculation and the class,
    call the perform function once, then get your desired results
    from the corresponding arrays.
    Access to the arrays is provided by get_ functions.

    \ingroup libwfa
 **/
class ex_analyse_ad {
private:

    double rh[3][2]; //!< \f$<r_h>\f$ for x,y,z and alpha/beta
    double rh2[3][2];//!< \f$<(r_h)²>\f$ for x,y,z and alpha/beta
    double re[3][2];//!< \f$<r_e>\f$ for x,y,z and alpha/beta
    double re2[3][3];//!< \f$<(r_e)²>\f$ for x,y,z and alpha/beta
    double sep[2];//!< \f$|<r_h - r_e>|\f$ for alpha/beta
    double sig_h[2];//!< \f$sqrt(<(r_h)²>-<r_h>²)\f$ for alpha/beta
    double sig_e[2];//!< \f$sqrt(<(r_e)²>-<r_e>²)\f$ for alpha/beta
    double prom[2];//!< Promotion number
    bool m_aeqb;//!< is alpha equal to beta

    /** \brief Forms all needed expected values over all spacial coordinates
     and spins.
     \param att Attachment matrix
     \param det Detachment matrix
     \param op Contract interface reference
     **/
    void ex_form_ad(const ab_matrix &att,
    		const ab_matrix &det, const multipol_con_i &op);
    /** \brief Calculates the separation \f$|<r_h - r_e>|\f$.
     \param spin Char to choose the spin --> a,b
     \return result
     **/
    double ex_mean_sep_ad(char spin);
    /** \brief Calculates \f$(\OMEGA)_h=sqrt(<(r_h)²>-<r_h>²)\f$.
     \param spin Char to choose the spin --> a,b
     \return result
     **/
    double ex_sig_h_ad (char spin);
    /** \brief Calculates the \f$(\OMEGA)_e=sqrt(<(r_e)²>-<r_e>²)\f$.
     \param spin Char to choose the spin --> a,b
     \return result
     **/
    double ex_sig_e_ad (char spin);

public:
    /** \brief Returns the exp. value of <r_h> for a sp. coordinate and spin.
     \param coord Char to choose the sp. coord.--> x,y,z
     \param spin Char to choose the spin --> a,b
     **/
    double get_rh (char coord, char spin);
    /** \brief Returns the exp. value of <(r_h)²> for a sp. coordinate and spin.
     \param coord Char to choose the sp. coord.--> x,y,z
     \param spin Char to choose the spin --> a,b
     **/
    double get_rh2 (char coord, char spin);
    /** \brief Returns the exp. value of <r_e> for a sp. coordinate and spin.
     \param coord Char to choose the sp. coord.--> x,y,z
     \param spin Char to choose the spin --> a,b
     **/
    double get_re (char coord, char spin);
    /** \brief Returns the exp. value of <(r_e)²> for a sp. coordinate and spin.
     \param coord Char to choose the sp. coord.--> x,y,z
     \param spin Char to choose the spin --> a,b
     **/
    double get_re2 (char coord, char spin);
    /** \brief Returns the separation for a sp. spin.
     \param spin Char to choose the spin --> a,b
     **/
    double get_sep (char spin);
    /** \brief Returns the sigma for the hole for a sp. spin.
     \param spin Char to choose the spin --> a,b
     **/
    double get_sig_h (char spin);
    /** \brief Returns the sigma for the electron for a sp. spin.
     \param spin Char to choose the spin --> a,b
     **/
    double get_sig_e (char spin);
    /** \brief Performs every analysis in this class and saves it in arrays.
     \param att Attachment matrix
     \param det Detachment matrix
     \param op Contract interface reference
     **/
    void perform (const ab_matrix &att, const ab_matrix &det,
    		const multipol_con_i &op);
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
};

}//end namespace libwfa

#endif /* EX_ANALYSE_AD_H_ */
