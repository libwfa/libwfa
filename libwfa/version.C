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


std::string version::get_license() {

    std::ostringstream ss;
    ss << "Copyright (c) 2014, F. Plasser and M. Wormit" << std::endl;
    ss << "All rights reserved." << std::endl << std::endl;
    ss << "Redistribution and use in source and binary forms, with or without"
            << std::endl;
    ss << "modification, are permitted provided that the following conditions"
            << " are met:" << std::endl << std::endl;
    ss << "1. Redistributions of source code must retain the above copyright"
            << " notice, this" << std::endl;
    ss << "   list of conditions and the following disclaimer." << std::endl;
    ss << "2. Redistributions in binary form must reproduce the above"
            << " copyright notice," << std::endl;
    ss << "   this list of conditions and the following disclaimer in the"
            << " documentation" << std::endl;
    ss << "   and/or other materials provided with the distribution."
            << std::endl << std::endl;
    ss << "3. Neither the name of the copyright holder nor the names of its"
            << " contributors" << std::endl;
    ss << "   may be used to endorse or promote products derived from this"
            << " software" << std::endl;
    ss << "   without specific prior written permission."
            << std::endl << std::endl;
    ss << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS"
            << " \"AS IS\" AND" << std::endl;
    ss << "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,"
            << " THE IMPLIED" << std::endl;
    ss << "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE"
            << " ARE" << std::endl;
    ss << "DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS"
            << " BE LIABLE" << std::endl;
    ss << "FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR"
            << " CONSEQUENTIAL" << std::endl;
    ss << "DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE"
            << " GOODS OR" << std::endl;
    ss << "SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)"
            << " HOWEVER" << std::endl;
    ss << "CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT"
            << " LIABILITY, " << std::endl;
    ss << "OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT"
            << " OF THE USE " << std::endl;
    ss << "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH"
            << " DAMAGE." << std::endl;

    return ss.str();
}

} // namespace libwfa

