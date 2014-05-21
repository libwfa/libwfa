#include <algorithm>
#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "pop_printer_default.h"

namespace libwfa {


void pop_printer_default::perform(const pop_data &p) {

    if (p.size() == 0) return;

    // Check if all data sets have the correct size
    for (pop_data::iterator i = p.begin(); i != p.end(); i++) {
        if (m_labels.size() != p.data(i).size()) {
            throw libwfa_exception("pop_printer_default",
                    "perform(const pop_data&)", __FILE__, __LINE__,
                    "Length of population data.");
        }
    }

    // Compute max width of an index number
    size_t nw1 = 0, nl = m_labels.size();
    while (nl != 0) { nl /= 10; nw1++; }

    // Compute max width of a label
    size_t nw2 = 0;
    for (std::vector<std::string>::const_iterator i = m_labels.begin();
            i != m_labels.end(); i++) {

        nw2 = std::max(nw2, i->size());
    }

    // Compute width of initial columns
    size_t nw = nw1 + nw2 + 1;
    nw = std::max(nw, (size_t) 4);

    std::string offset(nw - (nw1 + nw2 + 1), ' ');

    // Reduce column width, if total width is wider than 80 characters
    size_t maxwidth = 80, colwidth = m_colwidth;
    size_t mincolwidth = m_prec + 4;
    size_t width = nw + colwidth * p.size();
    while (width > maxwidth && colwidth > mincolwidth) {
        width -= p.size();
        colwidth--;
    }

    // Print header
    m_out << std::setw(nw) << std::right << "Atom";
    for (pop_data::iterator i = p.begin(); i != p.end(); i++) {
        m_out << std::setw(colwidth) << std::right << p.name(i);
    }
    m_out << std::endl;
    m_out << std::string(width, '-') << std::endl;
        
    std::vector<double> total(p.size(), 0.0);
    for (size_t i = 0, j = 1; i != m_labels.size(); i++, j++) {

        m_out << offset << std::setw(nw1) << j;
        m_out << " " << std::setw(nw2) << m_labels[i];

        size_t k = 0;
        for (pop_data::iterator kk = p.begin();
                kk != p.end(); k++, kk++) {

            const std::vector<double> &set = p.data(kk);
            m_out << std::setw(colwidth) << std::right << std::fixed <<
                    std::setprecision(m_prec) << set[i];
            total[k] += set[i];
        }
        m_out << std::endl;
    }

    // sum
    m_out << std::string(width, '-') << std::endl;

    m_out << std::setw(nw) << std::right << "Sum:";
    for (std::vector<double>::const_iterator i = total.begin();
            i != total.end(); i++) {
        m_out << std::setw(colwidth) << std::right << *i;
    }
    m_out << std::endl;
}


} // namespace libwfa

