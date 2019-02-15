#include "ctnum_analysis.h"
#include "om_descriptor.h"

namespace libwfa {

using namespace arma;


ctnum_analysis::ctnum_analysis(const mat &s, const uvec &b2p, const std::string &method,
                               const std::vector<std::string> &prop_list,
                               const std::vector<std::vector<int>> &at_lists):
                               m_nparts(0), m_s(s), m_b2p(b2p), ctnum_method(method),
                               prop_list(prop_list), at_lists(at_lists) {

    m_nparts = b2p.max() + 1;

    if (ctnum_method == "lowdin") {
        m_s_sqrt = sqrtmat_sympd(m_s);
    }
}


void ctnum_analysis::perform(const mat &tdm, mat &om) const {

    mat om_ao;
    form_om(m_s, m_s_sqrt, tdm, ctnum_method, om_ao);

    om.resize(m_nparts, m_nparts);
    om.fill(0.0);

    for (size_t i = 0; i < m_b2p.size(); i++) {

        size_t iat = m_b2p[i];
        for (size_t j = 0; j < m_b2p.size(); j++) {

            size_t jat = m_b2p[j];
            om(iat, jat) += om_ao(i, j);
        }
    }
}


void ctnum_analysis::form_om(const arma::mat &s, const arma::mat &s_sqrt,
    const arma::mat &tdm, const std::string &method, arma::mat &om) {

    if (method == "mulliken") {
        om = 0.5 * ((tdm * s) % (s * tdm) + tdm % (s * tdm * s));
    }
    else if (method == "lowdin") {
        om = square(s_sqrt * tdm * s_sqrt);
    }

}


auto ctnum_analysis::compute_omAt(const arma::mat &om, const std::vector<std::vector<int>> &blocks,
                                        const int &natoms) {
    mat om_at(natoms, natoms);
    om_at.fill(0.0);

    for (auto const& iblock: blocks) {
        for (auto const& jblock: blocks) {
            om_at(iblock[0], jblock[0]) = accu( om(span(iblock[1], iblock[2]), span(jblock[1], jblock[2])) );
        }
    }

    return om_at;
}


auto ctnum_analysis::bf_blocks(int &num_bas) {

    std::vector<std::vector<int>> blocks;

    return blocks;
}


auto ctnum_analysis::compute_omFrag(const arma::mat &om_at, const std::vector<std::vector<int>> &at_lists) {

    mat om_frag(at_lists.size(), at_lists.size());
    om_frag.fill(0.0);

    for (size_t i = 0; i < at_lists.size(); i++) {
        for (size_t j = 0; j < at_lists.size(); j++) {
            for (auto const& iatom: at_lists[i]) {
                for (auto const& jatom: at_lists[j]) {

                    om_frag(i, j) += om_at(iatom-1, jatom-1);

                }
            }
        }
    }

    return om_frag;

}


auto ctnum_analysis::compute_desc(const double &om_tot, const mat &om_frag, std::vector<std::string> &prop_list) {

    OmDescriptor desc(om_tot, om_frag);
    desc.ret_desc(prop_list);
    return desc.descriptor;

}


} // namespace libwfa
