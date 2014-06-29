#ifndef LIBWFA_VERSION_H
#define LIBWFA_VERSION_H

#include <list>
#include <string>

namespace libwfa {


/** \brief Version of the wave function analysis library

    The %version number consists of
    - major number
    - minor number
    - status string describing the release status

    For example, %version 2.0-alpha2 has major number 2, minor number 0,
    and status "alpha2" meaning the second alpha release.

    TODO:
    - add patch number
    - use cmakes version determination functions to set the version numbers

    \ingroup libwfa
 **/
class version {
private:
    static const unsigned k_major = 1; //!< Major %version number
    static const unsigned k_minor = 0; //!< Minor %version number
    static const char *k_status; //!< Version status
    static const char *k_authors[]; //!< List of authors

public:
    /** \brief Returns the major %version number
     **/
    static unsigned get_major() {
        return k_major;
    }

    /** \brief Returns the minor %version number
     **/
    static unsigned get_minor() {
        return k_minor;
    }

    /** \brief Returns the %version status
     **/
    static std::string get_status() {
        return std::string(k_status);
    }

    /** \brief Returns the string that corresponds to the %version
     **/
    static std::string get_string();

    /** \brief Returns a list of authors
     **/
    static std::list<std::string> get_authors();
};


} // namespace libwfa

#endif // LIBWFA_VERSION_H
