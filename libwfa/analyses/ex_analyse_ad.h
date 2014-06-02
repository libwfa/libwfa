#ifndef EX_ANALYSE_AD_H_
#define EX_ANALYSE_AD_H_

#include <libwfa/core/ab_matrix.h>

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
    /** \brief Calculates the expected value for the multipol of an
     exciton and returns either the alpha or beta value.
     \param[in] ad Either attachment or detachment matrix
     \param[in] op Multipoloperator
     \param[in] om Omega Matrix
     \param[in] spin Char Switch to set either 'a'=alpha or 'b'=beta
     **/
    double ex_multip_ad(const ab_matrix &ad, const ab_matrix &op,
            const ab_matrix &om, char spin);

    /** \brief Forms all needed expected values over all spacial coordinates
     and spins.
     \param[in] att Attachment matrix
     \param[in] det Detachment matrix
     \param[in] mxx Multipoloperator for r² in x coordinate
     \param[in] mx Multipoloperator for r in x coordinate
     \param[in] myy Multipoloperator for r² in y coordinate
     \param[in] my Multipoloperator for r in y coordinate
     \param[in] mzz Multipoloperator for r² in z coordinate
     \param[in] mz Multipoloperator for r in z coordinate
     \param[in] om Omega Matrix
     \param[out] rh[3][2] Multidimensional array to save
     ex. values of <rh> --> 3 coordinates 2 spins.
     \param[out] re[3][2] Multidimensional array to save
     ex. values of <re> --> 3 coordinates 2 spins.
     \param[out] rh2[3][2] Multidimensional array to save
     ex. values of <rh²> --> 3 coordinates 2 spins.
     \param[out] re2[3][2] Multidimensional array to save
     ex. values of <re²> --> 3 coordinates 2 spins.
     **/
    double ex_form_ad(const ab_matrix &att, const ab_matrix &det, const ab_matrix &mx,
            const ab_matrix &mxx, const ab_matrix &my, const ab_matrix &myy,
            const ab_matrix &mz, const ab_matrix &mzz, const ab_matrix &om);
    /** \brief Calculates the <rh>-<re> vector and returns its quantity.
     \param[in] spin Char to choose the spin --> a,b
     **/
    double ex_mean_sep_ad(char spin);
    /** \brief Calculates the sigma for the hole.
     \param[in] spin Char to choose the spin --> a,b
     **/
    double ex_sig_h_ad (char spin);
    /** \brief Calculates the sigma for the electron.
     \param[in] spin Char to choose the spin --> a,b
     **/
    double ex_sig_e_ad (char spin);

public:
    /** \brief Returns the exp. value of <rh> for a sp. coordinate and spin.
     \param[in] coord Char to choose the sp. coord.--> x,y,z
     \param[in] spin Char to choose the spin --> a,b
     **/
    double get_rh (char coord, char spin);
    /** \brief Returns the exp. value of <rh2> for a sp. coordinate and spin.
     \param[in] coord Char to choose the sp. coord.--> x,y,z
     \param[in] spin Char to choose the spin --> a,b
     **/
    double get_rh2 (char coord, char spin);
    /** \brief Returns the exp. value of <re> for a sp. coordinate and spin.
     \param[in] coord Char to choose the sp. coord.--> x,y,z
     \param[in] spin Char to choose the spin --> a,b
     **/
    double get_re (char coord, char spin);
    /** \brief Returns the exp. value of <re2> for a sp. coordinate and spin.
     \param[in] coord Char to choose the sp. coord.--> x,y,z
     \param[in] spin Char to choose the spin --> a,b
     **/
    double get_re2 (char coord, char spin);
    /** \brief Returns the separation for a sp. spin.
     \param[in] spin Char to choose the spin --> a,b
     **/
    double get_sep (char spin);
    /** \brief Returns the sigma for the hole for a sp. spin.
     \param[in] spin Char to choose the spin --> a,b
     **/
    double get_sig_h (char spin);
    /** \brief Returns the sigma for the electron for a sp. spin.
     \param[in] spin Char to choose the spin --> a,b
     **/
    double get_sig_e (char spin);
    /** \brief Performs every analysis in this class and saves it in arrays.
     \param[in] att Attachment matrix
     \param[in] det Detachment matrix
     \param[in] mxx Multipoloperator for r² in x coordinate
     \param[in] mx Multipoloperator for r in x coordinate
     \param[in] myy Multipoloperator for r² in y coordinate
     \param[in] my Multipoloperator for r in y coordinate
     \param[in] mzz Multipoloperator for r² in z coordinate
     \param[in] mz Multipoloperator for r in z coordinate
     \param[in] om Omega Matrix
    \param[out] sep Electron hole seperation
    \param[out] sig_h Varianz (sigma) for the hole
    \param[out] sig_e Varianz (sigma) for the electron
     **/
    void perform (const ab_matrix &att, const ab_matrix &det,
            const ab_matrix &mx,const ab_matrix &mxx, const ab_matrix &my,
            const ab_matrix &myy, const ab_matrix &mz, const ab_matrix &mzz,
            const ab_matrix &om);

};

}//end namespace libwfa

#endif /* EX_ANALYSE_AD_H_ */
