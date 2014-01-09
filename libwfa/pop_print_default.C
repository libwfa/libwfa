#include <algorithm>
#include <iomanip>
#include "pop_print_default.h"

namespace libwfa {


void pop_print_default::perform(const pop_data &p) {

    if (p.size() == 0) return;

    // Check if all data sets have the correct size
    for (pop_data::const_iterator i = p.begin(); i != p.end(); i++) {
        if (m_labels.size() != i->second.size()) {
            throw ;
        }
    }

    // Compute max width of an index number
    size_t nw1 = 0, nl = m_labels.size();
    while (nl != 0) { nl %= 10; nw1++; }

    // Compute max width of a label
    size_t nw2 = 0;
    for (std::vector<std::string>::const_iterator i = m_labels.begin();
            i != m_labels.end(); i++) {

        nw2 = std::max(nw2, i->size());
    }

    // Reduce column width, if total width is wider than 80 characters
    size_t maxwidth = 80, colwidth = m_colwidth;
    size_t mincolwidth = m_prec + 4;
    size_t width = nw1 + nw2 + 1 + colwidth * p.size();
    while (width > maxwidth && colwidth > mincolwidth) {
        width -= p.size();
        colwidth--;
    }

    // Print header
    m_out << std::setw(nw1 + nw2 + 1) << std::right << "Atom";
    for (pop_data::const_iterator i = p.begin(); i != p.end(); i++) {
        m_out << std::setw(colwidth) << std::right << i->first;
    }
    m_out << std::endl;
    m_out << std::string(width, '-') << std::endl;
        
    std::vector<double> total(p.size(), 0.0);
    for (size_t i = 0, j = 1; i != m_labels.size(); i++, j++) {

        m_out << std::setw(nw1) << j << " ";
        m_out << std::setw(nw2) << m_labels[i];

        size_t k = 0;
        for (pop_data::const_iterator kk = p.begin();
                kk != p.end(); k++, kk++) {

            m_out << std::setw(colwidth) << std::setprecision(m_prec) <<
                    std::right << kk->second[i];
            total[k] += kk->second[i];
        }
        m_out << std::endl;
    }

    // sum
    m_out << std::string(width, '-') << std::endl;

    m_out << std::setw(nw1 + nw2 + 1) << "Sum:";
    for (std::vector<double>::const_iterator i = total.begin();
            i != total.end(); i++) {
        m_out << std::setw(colwidth) << std::right << *i;
    }
    m_out << std::endl;
}


} // namespace libwfa

