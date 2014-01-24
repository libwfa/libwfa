#ifndef LIBWFA_NTO_DATA_H
#define LIBWFA_NTO_DATA_H

#include <iostream>
#include <list>
#include "ab_vector.h"
#include "dm_type.h"

namespace libwfa {

/** \brief Interface class to extract information on occupation numbers the
        eigenvector of a density matrix

    \ingroup libwfa
 **/
class nto_data_i {
public:
    virtual ~nto_data_i() { }

    /** \brief Perform the operation
        \param type Type of density matrix
        \param ni Occupation number vector
        \return Number of important occupation numbers
     **/
    virtual size_t perform(dm_type type, const ab_vector &ni) = 0;
};


/** \brief Implementation of nto_data_i to print to output stream

    \ingroup libwfa
 **/
class nto_data_print : public nto_data_i {
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
    nto_data_print(std::ostream &out, double thresh, size_t nnto) :
        m_out(out), m_thresh(thresh), m_nnto(nnto) { }

    /** \copydoc nto_data_i::perform
     **/
    virtual size_t perform(dm_type type, const ab_vector &ni);
};


/** \brief Implementation of nto_data_i to store the data

    \ingroup libwfa
 **/
class nto_data_extract : public nto_data_i {
public:
    struct ntoinfo {
        dm_type type;
        double total;
        double pr;
        std::vector<double> ni;

        ntoinfo(size_t n = 0) :
            type(dm_type::sdm), total(0.0), pr(0.0), ni(n, 0.0) { }
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
