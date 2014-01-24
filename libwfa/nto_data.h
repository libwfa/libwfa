#ifndef LIBWFA_NTO_DATA_H
#define LIBWFA_NTO_DATA_H

#include <iostream>
#include <list>
#include "ev_data_i.h"

namespace libwfa {

/** \brief Implementation of nto_data_i to print to output stream

    \ingroup libwfa
 **/
class nto_data_print : public ev_data_i {
public:
    static const char k_clazz[]; //!< Class name

private:
    std::ostream &m_out; //!< Output stream
    double m_thresh; //!< Threshold of important NTOs
    size_t m_nnto; //!< Number of leading occupation numbers to print

public:
    /** \brief Constructor
        \param out Output stream
        \param thresh Threshold for important NTOs
        \param nnto # of leading occupation numbers to print
     */
    nto_data_print(std::ostream &out, double thresh = 1e-6, size_t nnto = 3) :
        m_out(out), m_thresh(thresh), m_nnto(nnto) { }

    /** \copydoc nto_data_i::perform
     **/
    virtual size_t perform(dm_type type, const ab_vector &ni);

private:
    size_t print(const arma::Col<double> &ni);
};


/** \brief Implementation of nto_data_i to store the data

    \ingroup libwfa
 **/
class nto_data_extract : public ev_data_i {
public:
    static const char k_clazz[]; //!< Class name

public:
    struct ntoinfo {
        dm_type type;
        double total;
        double pr;
        std::vector<double> ni;

        ntoinfo(size_t n = 0) :
            type(dm_type::state), total(0.0), pr(0.0), ni(n, 0.0) { }
    };
    typedef ab_object<ntoinfo> ab_ntoinfo;
    typedef std::list<ab_ntoinfo>::const_iterator iterator;

private:
    double m_thresh; //!< Threshold of important NTOs
    size_t m_nnto; //!< Number of leading occupation numbers to print
    std::list<ab_ntoinfo> m_sets;

public:
    nto_data_extract(double thresh, size_t nnto) :
        m_thresh(thresh), m_nnto(nnto) { }

    virtual size_t perform(dm_type type, const ab_vector &ni);

    void clear() { return m_sets.clear(); }
    size_t size() const { return m_sets.size(); }
    iterator begin() const { return m_sets.begin(); }
    iterator end() const { return m_sets.end(); }
    const ab_ntoinfo &get_data(iterator i) { return *i; }
};


} // namespace adcman

#endif
