#include "libwfa_exception.h"
#include "density_list.h"

namespace libwfa {


const density_type density_type::state_density = density_type(0);
const density_type density_type::trans_density = density_type(1);
const density_type density_type::att_density = density_type(2);
const density_type density_type::det_density = density_type(3);


density_list::~density_list() {
    while (! m_lst.empty()) {
        delete m_lst.back();
        m_lst.pop_back();
    }
}


void density_list::add(size_t idx, const ab_matrix &dm) {

    std::list<dm_container *>::iterator i = m_lst.begin();
    for (; i != m_lst.end(); i++) if ((*i)->idx == idx) break;

    if (i != m_lst.end()) {
        throw libwfa_exception("density_list",
                "add(const ab_matrix &, size_t)", __FILE__, __LINE__,
                "Duplicate index.");
    }

    m_lst.push_back(new dm_container(idx, dm));
}


} // namespace libwfa

