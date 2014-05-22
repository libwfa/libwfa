#ifndef LIBWFA_DENSITY_TYPE_H
#define LIBWFA_DENSITY_TYPE_H

#include <ostream>
#include <string>

namespace libwfa {

/** \brief Enum class of density matrix types

    \ingroup libwfa
 **/
class density_type {
private:
    char m_type;

private:
    explicit density_type(char t = 0) : m_type(t) { }

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
};


std::ostream &operator<<(std::ostream &out, density_type type);


} // namespace libwfa

#endif // LIBWFA_DENSITY_TYPE_H
