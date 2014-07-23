#include <sstream>
#include "version.h"

namespace libwfa {

const char *version::k_status = "trunk";
const char *version::k_authors[] = {
        "Felix Plasser",
        "Stefanie Baeppler",
        "Benjamin Thomitzni",
        "Michael Wormit"
};


std::string version::get_string() {

    std::ostringstream ss;
    ss << k_major << "." << k_minor << "-" << std::string(k_status);
    return ss.str();
}


std::list<std::string> version::get_authors() {

    std::list<std::string> authors;
    size_t nauthors = sizeof(k_authors) / sizeof(char*);
    for(size_t i = 0; i < nauthors; i++) {
        authors.push_back(std::string(k_authors[i]));
    }

    return authors;
}


} // namespace libwfa

