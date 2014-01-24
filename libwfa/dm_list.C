#include "dm_list.h"
#include "libwfa_exception.h"

namespace libwfa {


const dm_type dm_type::sdm = dm_type(0);
const dm_type dm_type::tdm = dm_type(1);
const dm_type dm_type::adm = dm_type(2);
const dm_type dm_type::ddm = dm_type(3);
const dm_type dm_type::edm = dm_type(4);
const dm_type dm_type::hdm = dm_type(5);


void dm_list::add(size_t idx, const ab_matrix &dm) {

    if (index_exists(idx)) {
        throw libwfa_exception("dm_list", "add(size_t, const ab_matrix &)",
                __FILE__, __LINE__, "Duplicate index.");
    }

    m_lst[idx] = new ab_matrix(dm);
}


void dm_list::erase(size_t idx) {

    std::map<size_t, ab_matrix *>::iterator i = m_lst.find(idx);
    if (i != m_lst.end()) { m_lst.erase(i); }
}


void dm_list::clear() {

    for (std::map<size_t, ab_matrix *>::iterator i = m_lst.begin(); i != m_lst.end(); i++) {
        delete i->second; i->second = 0;
    }
    m_lst.clear();
}


} // namespace libwfa

