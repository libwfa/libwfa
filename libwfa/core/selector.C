#include "selector.h"


namespace libwfa {


void selector::select(size_t i) {

    if (i >= m_ntotal) return;
#ifdef LIBWFA_DEBUG
    if (m_indexes.count(i) != 0) {
        throw libwfa_exception("selector", "select(size_t)",
                __FILE__, __LINE__, "Already selected.");
    }
#endif
    m_indexes[i] = m_cur++;
}


void selector::select(size_t i0, size_t i1, size_t inc, bool reverse) {

    if (reverse)
        for (size_t i = i0, j = i1 - 1; i < i1; i += inc, j -= inc) select(j);
    else
        for (size_t i = i0; i < i1; i += inc) select(i);
}


void selector::select_all() {

    if (m_indexes.size() == 0) {
        m_cur = 0;
        for (size_t i = 0; i < m_ntotal; i++) m_indexes[i] = m_cur++;
    }
    else {
        std::map<size_t, size_t>::iterator it = m_indexes.begin();
        for (size_t i = 0; i < m_ntotal; i++) {
            if (i < it->first) m_indexes[i] = m_cur++;
            else it++;
        }
    }
}


void selector::deselect(size_t i) {

    if (i >= m_ntotal) return;

    std::map<size_t, size_t>::iterator it = m_indexes.find(i);
    if (it == m_indexes.end()) return;

    size_t pos = it->second;
    m_indexes.erase(it);
    for (it = m_indexes.begin(); it != m_indexes.end(); it++) {
        if (it->second > pos) it->second--;
    }
    m_cur--;
}


void selector::deselect_all() {

    m_indexes.clear();
    m_cur = 0;
}


std::vector<size_t> selector::get_selected() const {

    std::vector<size_t> el(m_indexes.size(), 0);
    for (std::map<size_t, size_t>::const_iterator i = m_indexes.begin();
            i != m_indexes.end(); i++) {
        el[i->second] = i->first;
    }
    return el;
}


arma::Col<arma::uword> selector::get_selected_arma() const {

    arma::Col<arma::uword> el(m_indexes.size());
    for (std::map<size_t, size_t>::const_iterator i = m_indexes.begin();
            i != m_indexes.end(); i++) {
        el(i->second) = i->first;
    }
    return el;
}


} // namespace libwfa

