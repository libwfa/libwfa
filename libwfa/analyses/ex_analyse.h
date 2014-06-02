#ifndef EX_ANALYSE_H_
#define EX_ANALYSE_H_

#include <libwfa/core/ab_matrix.h>

namespace libwfa {

/** \brief Class for methods to analyse excitons TDM.

    This class contains several methods to analyse exciton TDMs.
    You can:
    - calculate the expected value for the multipol operators and saves them
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
    /** \brief Calculates the expected value for the multipol of an
             exciton and returns either the alpha or beta value.
             \param[in] tdm Transition Density Matrix
             \param[in] op1 Multipoloperator for hole
             \param[in] op2 Multipoloperator for electron
             \param[in] om Omega Matrix
             \param[in] spin Char Switch to set either 'a'=alpha or 'b'=beta
             **/
    double ex_multip(const ab_matrix &tdm, const ab_matrix &op1,
            const ab_matrix &op2, const ab_matrix &om, char spin);
        /** \brief Forms all needed expected values over all spacial coordinates
         and spins.
         \param[in] tdm Transition Density Matrix
         \param[in] s Overlap matrix
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
         \param[out] rhre[2] Multidimensional array to save
         ex. values of <rhre> --> 3 coordinates 2 spins.
         **/
    void ex_form(const ab_matrix &tdm, const ab_matrix &s, const ab_matrix &mxx,
            const ab_matrix &mx, const ab_matrix &myy, const ab_matrix &my,
            const ab_matrix &mzz, const ab_matrix &mz, const ab_matrix &om);
        /** \brief Calculates the expected distance over one
         spacial coordinate and for a specific spin.
         \param[in] koord Char to choose the sp. coord.--> x,y,z
         \param[in] spin Char to choose the spin --> a,b
         **/
    double ex_d_ex_c(char koord, char spin);
        /** \brief Calculates the expected distance over all
         spacial coordinate and for a specific spin.
         \param[in] koord Char to choose the sp. coord.--> x,y,z
         \param[in] spin Char to choose the spin --> a,b
         **/
    double ex_d_ex_tot(char spin);
    /** \brief Calculates the <rh>-<re> vector and returns its quantity.
         \param[in] spin Char to choose the spin --> a,b
         **/
    double ex_mean_sep(char spin);
        /** \brief Calculates the sigma for the hole.
         \param[in] spin Char to choose the spin --> a,b
         **/
    double ex_sig_h(char spin);
        /** \brief Calculates the sigma for the electron.
         \param[in] spin Char to choose the spin --> a,b
         **/
    double ex_sig_e(char spin);
        /** \brief Calculates the covalenz of rh, re.
         \param[in] spin Char to choose the spin --> a,b
         **/
    double ex_cov(char spin);
        /** \brief Calculates the correlation of rh, re.
         \param[in] spin Char to choose the spin --> a,b
         **/
    double ex_corr(char spin);

    public:
    /** \brief Returns the exp. value of <rh> for a sp. coordinate and spin.
             \param[in] koord Char to choose the sp. coord.--> x,y,z
             \param[in] spin Char to choose the spin --> a,b
             **/
        double get_rh(char koord, char spin);
            /** \brief Returns the exp. value of <re> for a sp. coordinate and spin.
             \param[in] koord Char to choose the sp. coord.--> x,y,z
             \param[in] spin Char to choose the spin --> a,b
             **/
        double get_re(char koord, char spin);
            /** \brief Returns the exp. value of <rh2> for a sp. coordinate and spin.
             \param[in] koord Char to choose the sp. coord.--> x,y,z
             \param[in] spin Char to choose the spin --> a,b
             **/
        double get_rh2(char koord, char spin);
            /** \brief Returns the exp. value of <re2> for a sp. coordinate and spin.
             \param[in] koord Char to choose the sp. coord.--> x,y,z
             \param[in] spin Char to choose the spin --> a,b
             **/
        double get_re2(char koord, char spin);
            /** \brief Returns the exp. value of <rhre> for a sp. coordinate and spin.
             \param[in] koord Char to choose the sp. coord.--> x,y,z
             \param[in] spin Char to choose the spin --> a,b
             **/
        double get_rhre(char koord, char spin);
        /** \brief Performs every analysis in this class and saves it in arrays.
         \param[in] tdm Transition Density Matrix
         \param[in] s Overlap matrix
         \param[in] mxx Multipoloperator for r² in x coordinate
         \param[in] mx Multipoloperator for r in x coordinate
         \param[in] myy Multipoloperator for r² in y coordinate
         \param[in] my Multipoloperator for r in y coordinate
         \param[in] mzz Multipoloperator for r² in z coordinate
         \param[in] mz Multipoloperator for r in z coordinate
         \param[in] om Omega Matrix
         \param[out] dex_c Expected distance for one coordinate
         \param[out] dex_tot Total expected distance
         \param[out] sep Electron-hole seperation
        \param[out] sig_h Varianz (sigma) for the hole
        \param[out] sig_e Varianz (sigma) for the electron
        \param[out] cov Covariance of rh,re
        \param[out] corr Correlation factor of rh, re
         **/
        void perform(const ab_matrix &tdm, const ab_matrix &s, const ab_matrix &mxx,
                const ab_matrix &mx, const ab_matrix &myy, const ab_matrix &my,
                const ab_matrix &mzz, const ab_matrix &mz, const ab_matrix &om);
        /** \brief Returns the separation for a sp. spin.
             \param[in] spin Char to choose the spin --> a,b
             **/
        double get_sep(char spin);
        /** \brief Returns the distance for a sp. coordinate and spin.
             \param[in] koord Char to choose the sp. coord.--> x,y,z
             \param[in] spin Char to choose the spin --> a,b
             **/
        double get_dex_c (char coord, char spin);
        /** \brief Returns the distance for a sp. spin.
            \param[in] spin Char to choose the spin --> a,b
              **/
        double get_dex_tot (char spin);
        /** \brief Returns the sigma for the hole for a sp. spin.
            \param[in] spin Char to choose the spin --> a,b
              **/
        double get_sig_h (char spin);
        /** \brief Returns the sigma for the electron for a sp. spin.
            \param[in] spin Char to choose the spin --> a,b
             **/
        double get_sig_e (char spin);
        /** \brief Returns the covariance for a sp. spin.
            \param[in] spin Char to choose the spin --> a,b
             **/
        double get_cov (char spin);
        /** \brief Returns the correlation for a sp. spin.
            \param[in] spin Char to choose the spin --> a,b
             **/
        double get_corr (char spin);
}; // end Class

} // end namespace libwfa


#endif
