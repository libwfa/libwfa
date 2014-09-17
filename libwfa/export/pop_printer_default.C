#include <algorithm>
#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "pop_printer_default.h"

namespace libwfa {


void pop_printer_default::perform(const pop_data &p, std::ostream &out) const {

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

    std::string os1(m_offset, ' '), os2(nw - (nw1 + nw2 + 1) + m_offset, ' ');

    // Reduce column width, if total width is wider than 80 characters
    size_t maxwidth = 80, colwidth = m_colwidth;
    size_t mincolwidth = m_prec + 4;
    size_t width = m_offset + nw + colwidth * p.size();
    while (width > maxwidth && colwidth > mincolwidth) {
        width -= p.size();
        colwidth--;
    }

    // Print header
    out << std::right << std::fixed << std::setprecision(m_prec);
    out << os1 << std::setw(nw) << "Atom";
    for (pop_data::iterator i = p.begin(); i != p.end(); i++) 
        out << std::setw(colwidth) << p.name(i);

    out << std::endl;
    out << os1 << std::string(width, '-') << std::endl;
        
    arma::vec total(p.size(), arma::fill::zeros);
    for (size_t i = 0, j = 1; i != m_labels.size(); i++, j++) {

        out << os2 << std::setw(nw1) << j;
        out << " " << std::setw(nw2) << m_labels[i];

        size_t k = 0;
        for (pop_data::iterator kk = p.begin();
                kk != p.end(); k++, kk++) {

            const arma::vec &set = p.data(kk);
            total(k) += set(i);

            out << std::setw(colwidth) << set(i);
        }
        out << std::endl;
    }

    // sum
    out << os1 << std::string(width, '-') << std::endl;

    out << std::right;
    out << os1 << std::setw(nw) << "Sum:";
    for (size_t i = 0; i < total.size(); i++) 
        out << std::setw(colwidth) << total(i);
    out << std::endl << std::endl;
}


} // namespace libwfa

