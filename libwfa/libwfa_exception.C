#include <cstdio>
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

    if(strlen(m_clazz) == 0) {
        if(strlen(m_method) == 0) {
            if(strlen(m_file) == 0)
                snprintf(m_what, 1024, "libwfa: %s", m_message);
            else
                snprintf(m_what, 1024, "libwfa, %s (%u): %s",
                        m_file, m_line, m_message);
        } else {
            if(strlen(m_file) == 0)
                snprintf(m_what, 1024, "libwfa::%s: %s", m_method, m_message);
            else
                snprintf(m_what, 1024, "libwfa::%s, %s (%u): %s",
                        m_method, m_file, m_line, m_message);
        }
    } else {
        if(strlen(m_method) == 0) {
            if(strlen(m_file) == 0)
                snprintf(m_what, 1024, "libwfa::%s: %s",
                        m_clazz, m_message);
            else
                snprintf(m_what, 1024, "libwfa::%s, %s (%u): %s",
                        m_clazz, m_file, m_line, m_message);
        } else {
            if(strlen(m_file) == 0)
                snprintf(m_what, 1024, "libwfa::%s::%s: %s",
                        m_clazz, m_method, m_message);
            else
                snprintf(m_what, 1024, "libwfa::%s::%s, %s (%u): %s",
                        m_clazz, m_method, m_file, m_line, m_message);
        }
    }
}


const char *libwfa_exception::what() const throw() {

    return m_what;
}

} // namespace libwfa

