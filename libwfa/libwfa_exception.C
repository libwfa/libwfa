#include <sstream>
#include <cstring>
#include "libwfa_exception.h"

namespace libwfa {


void libwfa_exception::init(const char *clazz, const char *method,
    const char *file, unsigned int line, const char *message) throw() {

    if(clazz == NULL) m_clazz[0] = '\0';
    else { strncpy(m_clazz, clazz, 128); m_clazz[127] = '\0'; }
    if(method == NULL) m_method[0] = '\0';
    else { strncpy(m_method, method, 128); m_method[127] = '\0'; }
    if(file == NULL) m_file[0] = '\0';
    else { strncpy(m_file, file, 128); m_file[127] = '\0'; }
    m_line = line;
    if(message == NULL) m_message[0] = '\0';
    else { strncpy(m_message, message, 256); m_message[255] = '\0'; }

    if(strlen(m_message) == 0) {
        strcpy(m_message, "<No error message>");
    }

    std::ostringstream oss;
    if(strlen(m_clazz) == 0) {
        if(strlen(m_method) == 0) {
            if(strlen(m_file) == 0)
                oss << "libwfa: " << m_message;
            else
                oss << "libwfa, " << m_file << " (" << m_line << "): "
                    << m_message;
        } else {
            if(strlen(m_file) == 0)
                oss << "libwfa::" << m_method << ": " << m_message;
            else
                oss << "libwfa::" << m_method << ", " << m_file << " ("
                    << m_line << "): " << m_message;
        }
    } else {
        if(strlen(m_method) == 0) {
            if(strlen(m_file) == 0)
                oss << "libwfa::" << m_clazz << ": " << m_message;
            else
                oss << "libwfa::" << m_clazz << ", " << m_file << " ("
                    << m_line << "): " << m_message;
        } else {
            if(strlen(m_file) == 0)
                oss << "libwfa::" << m_clazz << "::" << m_method << ": "
                    << m_message;
            else
                oss << "libwfa::" << m_clazz << "::" << m_method << ", "
                    << m_file << " (" << m_line << "): " << m_message;
        }
    }
    strncpy(m_what, oss.str().c_str(), 1024);
}


const char *libwfa_exception::what() const throw() {

    return m_what;
}

} // namespace libwfa

