#ifndef LIBWFA_STATE_INFO_H
#define LIBWFA_STATE_INFO_H

namespace libwfa {

/** \brief Enum class for spin multiplicity

    \ingroup libwfa
 **/
class spin {
private:
    unsigned char m_s2; //!< Spin state

    spin(unsigned char s2) : m_s2(s2) { }

public:
    unsigned char multiplicity() const { return m_s2; }

    //! \name Comparison operators
    //@{
    bool operator==(const spin &s) const { return m_s2 == s.m_s2; }
    bool operator!=(const spin &s) const { return m_s2 == s.m_s2; }
    bool operator< (const spin &s) const { return m_s2 <  s.m_s2; }
    bool operator<=(const spin &s) const { return m_s2 <= s.m_s2; }
    bool operator> (const spin &s) const { return m_s2 >  s.m_s2; }
    bool operator>=(const spin &s) const { return m_s2 >= s.m_s2; }
    //@}

    //! \name Spin multiplicity
    //@{
    static const spin unspecified; //!< Spin of unrestricted calculation
    static const spin singlet; //!< Singlet spin
    static const spin doublet; //!< Doublet spin
    static const spin triplet; //!< Triplet spin
    static const spin quartet; //!< Quartet spin
    static const spin quintet; //!< Quintet spin
    //@}
};

/** \brief Information on the state

    \ingroup libwfa
 **/
struct state_info {
    unsigned int stateno; //!< Number of state within current irrep
    spin multiplicity; //!< Multiplicity
    std::string irrep; //!< Irreducible representation of state (not transition)
    double energy; //!< Excitation energy (in eV)
    double osc_strength; //!< Oscillator strength

    state_info(unsigned int n, spin m,
        const std::string &ir, double e, double os) :
        stateno(n), multiplicity(m), irrep(ir), energy(e), osc_strength(os) { }


    /** \brief Convert state information into term symbol
        \param info State information
        \param sep Separation character (default: space)
        \return String containing the term symbol

        Creates a string with the term symbol in the format:
        - "1_3_Bu" for first triplet state of irreducible representation Bu
        - "2_A1" for second state of irreducible representation A1 (if unrestricted)
        if '_' was chosen as separator.

     **/
    std::string convert(char sep = ' ') const;
};


} // namespace libwfa

#endif // LIBWFA_STATE_INFO_H
