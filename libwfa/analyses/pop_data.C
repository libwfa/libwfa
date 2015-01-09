#include <algorithm>
#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "pop_data.h"

namespace libwfa {


void pop_data::print(std::ostream &out, const std::vector<std::string> &l,
        size_t cw, size_t prec, size_t off) const {

    if (m_sets.size() == 0) return;

    // Check if all data sets have the correct size
    for (iterator i = m_sets.begin(); i != m_sets.end(); i++) {
        if (l.size() != i->data.size()) {
            throw libwfa_exception("pop_data", "print(...) const",
                    __FILE__, __LINE__, "Length of population data.");
        }
    }

    // Compute max width of an index number
    size_t nw1 = 0, nl = l.size();
    while (nl != 0) { nl /= 10; nw1++; }

    // Compute max width of a label
    size_t nw2 = 0;
    for (std::vector<std::string>::const_iterator i = l.begin();
            i != l.end(); i++) {
        nw2 = std::max(nw2, i->size());
    }

    // Compute width of initial columns
    size_t nw = nw1 + nw2 + 1;
    nw = std::max(nw, (size_t) 4);

    std::string os1(off, ' '), os2(nw - (nw1 + nw2 + 1) + off, ' ');

    // Reduce column width, if total width is wider than 80 characters
    size_t maxwidth = 80;
    size_t mincolwidth = prec + 4;
    size_t width = off + nw + cw * m_sets.size();
    while (width > maxwidth && cw > mincolwidth) {
        width -= m_sets.size();
        cw--;
    }

    // Print header
    out << std::right << std::fixed << std::setprecision(prec);
    out << os1 << std::setw(nw) << "Atom";
    for (pop_data::iterator i = m_sets.begin(); i != m_sets.end(); i++) {
        out << std::setw(cw) << i->name;
    }
    out << std::endl;

    out << os1 << std::string(width - off, '-') << std::endl;
        
    arma::vec total(m_sets.size(), arma::fill::zeros);
    for (size_t i = 0, j = 1; i != l.size(); i++, j++) {

        out << os2 << std::setw(nw1) << j;
        out << " " << std::setw(nw2) << l[i];

        size_t k = 0;
        for (iterator kk = m_sets.begin(); kk != m_sets.end(); k++, kk++) {

            total(k) += kk->data(i);

            out << std::setw(cw) << kk->data(i);
        }
        out << std::endl;
    }

    // sum
    out << os1 << std::string(width - off, '-') << std::endl;

    out << std::right;
    out << os1 << std::setw(nw) << "Sum:";
    for (size_t i = 0; i < total.size(); i++) 
        out << std::setw(cw) << total(i);
    out << std::endl;
}


} // namespace libwfa

