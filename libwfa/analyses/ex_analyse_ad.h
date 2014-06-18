#ifndef EX_ANALYSE_AD_H_
#define EX_ANALYSE_AD_H_

#include <libwfa/core/ab_matrix.h>
#include <libwfa/core/contract_ad_i.h>
#include <libwfa/core/contract_ad.h>

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

    /** \brief Forms all needed expected values over all spacial coordinates
     and spins.
     \param[in] att Attachment matrix
     \param[in] det Detachment matrix
     \param[in] s Overlap matrix
     \param[in] name Contract interface reference

     **/
    void ex_form_ad(const ab_matrix &att, const ab_matrix &det,
            const Mat<double> &s, const contract_ad_i &name);
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
     \param[in] s Overlap matrix
     \param[in] name Contract Interface reference
     **/
    void perform (const ab_matrix &att, const ab_matrix &det,
            const Mat<double> &s, const contract_ad_i &name);

};

}//end namespace libwfa

#endif /* EX_ANALYSE_AD_H_ */
