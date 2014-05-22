#include <libwfa/libwfa_exception.h>
#include "selector.h"


namespace libwfa {


void selector::select(size_t i) {

    check(i);
    if (! m_indexes[i]) m_nselected++;
    m_indexes[i] = true;
}


void selector::select(size_t i0, size_t i1, size_t inc) {

    check(i0);
    check(i1);
    for (size_t i = i0; i <= i1; i += inc) {
        if (! m_indexes[i]) m_nselected++;
        m_indexes[i] = true;
    }
}


void selector::select_all() {

    m_indexes = std::vector<bool>(m_indexes.size(), true);
    m_nselected = m_indexes.size();
}


void selector::deselect(size_t i) {

    check(i);
    if (m_indexes[i]) m_nselected--;
    m_indexes[i] = false;
}


void selector::deselect_all() {

    m_indexes = std::vector<bool>(m_indexes.size(), false);
    m_nselected = 0;
}


std::vector<size_t> selector::get_selected() const {

    std::vector<size_t> el(m_nselected, 0);
    size_t pos = 0;
    for (size_t i = 0; i < m_indexes.size() && pos < el.size(); i++) {
        if (m_indexes[i]) el[pos++] = i;
    }
    return el;
}


arma::Col<arma::uword> selector::get_selected_arma() const {

    arma::Col<arma::uword> el(m_nselected);
    size_t pos = 0;
    for (size_t i = 0; i < m_indexes.size() && pos < el.size(); i++) {
        if (m_indexes[i]) el[pos++] = i;
    }
    return el;
}


void selector::check(size_t i) const {

    if (i >= m_indexes.size()) {
        throw libwfa_exception("selector", "check(size_t &) const",
                __FILE__, __LINE__, "i");
    }
}



} // namespace libwfa

