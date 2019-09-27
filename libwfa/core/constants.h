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
    static const double au2rcm; //!< Hartree to cm^-1 factor
    static const double au2D;  //!< a.u. to Debye
    static const double alpha;  //!< fine-structure constant; 1/c0 in a.u.
    //@}
};

} // namespace libwfa

#endif // LIBWFA_CONVERSION_H
