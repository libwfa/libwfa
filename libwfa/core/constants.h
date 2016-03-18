#ifndef LIBWFA_CONVERSION_H
#define LIBWFA_CONVERSION_H

namespace libwfa {

/** \brief Constants and unit constants factors

    \ingroup base
 **/
struct constants {

    //! \name Length constants factors
    //@{
    static const double au2ang; //!< a.u. to Angstrom constants factor
    static const double au2nm; //!< a.u. to nm constants factor
    static const double au2eV; //!< Hartree to eV factor
    //@}
};

} // namespace libwfa

#endif // LIBWFA_CONVERSION_H
