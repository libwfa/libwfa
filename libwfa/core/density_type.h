#ifndef LIBWFA_DENSITY_TYPE_H
#define LIBWFA_DENSITY_TYPE_H

#include <bitset>
#include <ostream>
#include <string>

namespace libwfa {

/** \brief Enum class of density matrix types

    \ingroup libwfa
 **/
class density_type {
public:
    enum { NT = 7 };
    /** \brief Numeric values for density types to be used with flag_t
     **/
    enum dt {
        NONE = 0,       //!< NONE
        STATE  = 1 << 0,//!< STATE
        TRANS  = 1 << 1,//!< TRANS
        DIFF   = 1 << 2,//!< DIFF
        ATTACH = 1 << 3,//!< ATTACH
        DETACH = 1 << 4,//!< DETACH
        PART   = 1 << 5,//!< PART
        HOLE   = 1 << 6,//!< HOLE
        AD     = 24,    //!< AD
        PH     = 96,    //!< PH
        ALL    = 127    //!< ALL
    };
    typedef std::bitset<NT> flag_t; //!< Flag for density types

private:
    size_t m_type;

private:
    explicit density_type(size_t t = 0) : m_type(t) { }

public:
    //! \name Comparison operators
    //@{
    bool operator==(const density_type &t) const { return m_type == t.m_type; }
    bool operator!=(const density_type &t) const { return m_type != t.m_type; }
    bool operator< (const density_type &t) const { return m_type <  t.m_type; }
    bool operator<=(const density_type &t) const { return m_type <= t.m_type; }
    bool operator> (const density_type &t) const { return m_type >  t.m_type; }
    bool operator>=(const density_type &t) const { return m_type >= t.m_type; }
    //@}

    //! \name Density matrix types
    //@{
    static const density_type state; //!< State density matrix type
    static const density_type transition; //!< Transition density matrix type
    static const density_type difference; //!< Difference density matrix type
    static const density_type attach; //!< Attachment density matrix type
    static const density_type detach; //!< Detachment density matrix type
    static const density_type particle; //!< Particle density matrix type
    static const density_type hole; //!< Hole density matrix type
    //@}


    /** \brief Convert density_type into string
     **/
    std::string convert() const;

    /** \brief Convert density_type into a longer string
     **/
    std::string convert2() const;

    /** \brief Return whether density_type corresponds to a symmetric matrix
     **/
    bool is_symm() const;

    /** \brief Multiply wit this factor during export as cube file
     **/
    double export_factor() const;

    /** \brief Conversion to integer
     **/
    operator size_t() const { return m_type; }
};


std::ostream &operator<<(std::ostream &out, density_type type);


} // namespace libwfa

#endif // LIBWFA_DENSITY_TYPE_H
