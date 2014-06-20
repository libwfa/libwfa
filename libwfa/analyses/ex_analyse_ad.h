#ifndef EX_ANALYSE_AD_H_
#define EX_ANALYSE_AD_H_

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/contract_i.h>

namespace libwfa {

/** \brief Class for methods to analyse excitons through attach/detach matrices.

    This class contains several methods to analyse exciton TDMs.
    You can:
    - calculate the expected value for the multipol operators and saves them
    - calculate the mean seperation
    - calculate the sigma for either the hole or electron

    To use this class you will have to form the needed matrices,
    create an object for the class, call the perform function once,
    then get your desired results from the corresponding arrays.
    Access to the arrays is provided by get_ functions.

    \ingroup libwfa
 **/
class ex_analyse_ad {
private:

    double rh[3][2], rh2[3][2], re[3][2], re2[3][3];
    double sep[2];
    double sig_h[2];
    double sig_e[2];
    double prom[2];
    bool m_aeqb;

    /** \brief Forms all needed expected values over all spacial coordinates
     and spins.
     \param att Attachment matrix
     \param det Detachment matrix
     \param s Overlap matrix
     \param op Contract interface reference

     **/
    void ex_form_ad(const ab_matrix &att,
    		const ab_matrix &det, const contract_i &op);
    /** \brief Calculates the <rh>-<re> vector and returns its quantity.
     \param spin Char to choose the spin --> a,b
     \return result
     **/
    double ex_mean_sep_ad(char spin);
    /** \brief Calculates the sigma for the hole.
     \param spin Char to choose the spin --> a,b
     \return result
     **/
    double ex_sig_h_ad (char spin);
    /** \brief Calculates the sigma for the electron.
     \param spin Char to choose the spin --> a,b
     \return result
     **/
    double ex_sig_e_ad (char spin);

public:
    /** \brief Returns the exp. value of <rh> for a sp. coordinate and spin.
     \param coord Char to choose the sp. coord.--> x,y,z
     \param spin Char to choose the spin --> a,b
     **/
    double get_rh (char coord, char spin);
    /** \brief Returns the exp. value of <rh2> for a sp. coordinate and spin.
     \param coord Char to choose the sp. coord.--> x,y,z
     \param spin Char to choose the spin --> a,b
     **/
    double get_rh2 (char coord, char spin);
    /** \brief Returns the exp. value of <re> for a sp. coordinate and spin.
     \param coord Char to choose the sp. coord.--> x,y,z
     \param spin Char to choose the spin --> a,b
     **/
    double get_re (char coord, char spin);
    /** \brief Returns the exp. value of <re2> for a sp. coordinate and spin.
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
     \param s Overlap matrix
     \param op Contract interface reference
     **/
    void perform (const ab_matrix &att, const ab_matrix &det,
    		const contract_i &op);
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
