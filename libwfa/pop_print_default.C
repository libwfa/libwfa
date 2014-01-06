#include "pop_print_default.h"

namespace libwfa {


pop_print_default::pop_print_default(const std::vector<std::string> &l,
        std::ostream &out, size_t colwidth = 20, size_t prec = 6) :
    m_labels(l), m_out(out), m_colwidth(colwidth), m_prec(prec) {

}


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
            i != m_labels.end(); i++) { nw2 = std::max(nw2, i->size()); }

    // Reduce column width, if total width is wider than 80 characters
    size_t maxwidth = 80, colwidth = m_colwidth;
    size_t mincolwidth = m_prec + 4;
    size_t width = nw1 + nw2 + 1 + colwidth * m_pop.size();
    while (width > maxwidth && colwidth > mincolwidth) {
        width -= m_pop.size();
        colwidth--;
    }

    // Print header
    out << std::setw(nw1 + nw2 + 1) << std::right << "Atom";
    for (pop_data::const_iterator i = p.begin(); i != p.end(); i++) {
        out << std::setw(colwidth) << std::right << i->first;
    }
    out << std::endl;
    out << std::string(width, '-') << std::endl;
        
    std::vector<double> total(p.size(), 0.0);
    for (size_t i = 0, j = 1; i != m_labels.size(); i++, j++) {

            out << std::setw(nw1) << j << " ";
            out << std::setw(nw2) << m_labels[i];

            size_t k = 0;
            for (pop_data::const_iterator kk = p.begin();
                    kk != p.end(); k++, kk++) {

                out << std::setw(colwidth) << std::setprecision(m_prec) <<
                        std::right << kk->second[i];
                total[k] += kk->second[i];
            }
            out << std::endl;
        }

        // sum
        out << std::string(width, '-') << std::endl;

        out << std::setw(nw1 + nw2 + 1) << "Sum:";
        for (std::vector<double>::const_iterator i = total.begin();
                i != total.end(); i++) {
            out << std::setw(colwidth) << std::right << total[i];
        }
        out << std::endl;
    }
};


} // namespace adcman

#endif // ADCMAN_MULLIKEN_POPULATION_BASE_H
