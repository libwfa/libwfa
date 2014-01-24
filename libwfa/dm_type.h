#ifndef LIBWFA_DM_TYPE_H
#define LIBWFA_DM_TYPE_H

#include <string>

namespace libwfa {

/** \brief Enum class of density matrix types

    \ingroup libwfa
 **/
class dm_type {
private:
    char m_type;

private:
    explicit dm_type(char t = 0) : m_type(t) { }

public:
    //! \name Comparison operators
    //@{
    bool operator==(const dm_type &t) const { return m_type == t.m_type; }
    bool operator!=(const dm_type &t) const { return m_type != t.m_type; }
    bool operator< (const dm_type &t) const { return m_type <  t.m_type; }
    bool operator<=(const dm_type &t) const { return m_type <= t.m_type; }
    bool operator> (const dm_type &t) const { return m_type >  t.m_type; }
    bool operator>=(const dm_type &t) const { return m_type >= t.m_type; }
    //@}

    //! \name Density matrix types
    //@{
    static const dm_type sdm; //!< State density matrix type
    static const dm_type tdm; //!< Transition density matrix type
    static const dm_type adm; //!< Attachment density matrix type
    static const dm_type ddm; //!< Detachment density matrix type
    static const dm_type hdm; //!< Hole density matrix type
    static const dm_type edm; //!< Electron density matrix type
    //@}


    /** \brief Convert dm_type into string
     **/
    std::string convert() const;
};


} // namespace libwfa

#endif // LIBWFA_DM_TYPE_H
