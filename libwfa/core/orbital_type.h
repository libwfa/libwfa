#ifndef LIBWFA_ORBITAL_TYPE_H
#define LIBWFA_ORBITAL_TYPE_H

#include <bitset>
#include <ostream>
#include <string>

namespace libwfa {

/** \brief Enum class of density matrix types

    \ingroup libwfa
 **/
class orbital_type {
public:
    enum { NT = 4, OT_ALL = 15 };

    typedef std::bitset<NT> flag_t; //!< Flag for density types

private:
    size_t m_type;

private:
    explicit orbital_type(size_t t = 0) : m_type(t) { }

public:
    //! \name Comparison operators
    //@{
    bool operator==(const orbital_type &t) const { return m_type == t.m_type; }
    bool operator!=(const orbital_type &t) const { return m_type != t.m_type; }
    bool operator< (const orbital_type &t) const { return m_type <  t.m_type; }
    bool operator<=(const orbital_type &t) const { return m_type <= t.m_type; }
    bool operator> (const orbital_type &t) const { return m_type >  t.m_type; }
    bool operator>=(const orbital_type &t) const { return m_type >= t.m_type; }
    //@}

    //! \name Density matrix types
    //@{
    static const orbital_type mo; //!< Molecular orbital type
    static const orbital_type no; //!< Natural orbital type
    static const orbital_type nto; //!< Natural transition orbital type
    static const orbital_type ndo; //!< Natural difference orbital type
    //@}


    /** \brief Convert orbital_type into string
     **/
    std::string convert() const;

    /** \brief Conversion to integer
     **/
    operator size_t() { return m_type; }
};


std::ostream &operator<<(std::ostream &out, orbital_type type);


} // namespace libwfa

#endif // LIBWFA_ORBITAL_TYPE_H
