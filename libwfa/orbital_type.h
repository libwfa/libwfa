#ifndef LIBWFA_ORBITAL_TYPE_H
#define LIBWFA_ORBITAL_TYPE_H

#include <string>

namespace libwfa {

/** \brief Enum class of density matrix types

    \ingroup libwfa
 **/
class orbital_type {
private:
    char m_type;

private:
    explicit orbital_type(char t = 0) : m_type(t) { }

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
};


} // namespace libwfa

#endif // LIBWFA_ORBITAL_TYPE_H
