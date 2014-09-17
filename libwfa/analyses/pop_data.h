#ifndef LIBWFA_POP_DATA_H
#define LIBWFA_POP_DATA_H

#include <armadillo>
#include <list>
#include <string>

namespace libwfa {

/** \brief Container for several named population data sets

    \ingroup libwfa
 **/
class pop_data {
private:
    /** \brief Structure to keep one single data set
     **/
    struct data_set{
        std::string name;
        arma::vec data;

        data_set(const std::string &name_,
            const arma::vec &data_ = arma::vec()) :
            name(name_), data(data_) {
        }
    };

public:
    typedef std::list<data_set>::const_iterator iterator; //!< Iterator to traverse the data

private:
    std::list<data_set> m_sets; //!< Collection of data sets

public:
    /** \brief Add a new empty set for population data
        \param name Name of data set
        \return Reference to the empty data set.
     **/
    arma::vec &add(const std::string &name) {
        m_sets.push_back(data_set(name));
        return m_sets.back().data;
    }

    /** \brief Add a new set for population data
        \param name Name of data set
        \param data The data to be added
     **/
    void add(const std::string &name, const arma::vec &data) {
        m_sets.push_back(data_set(name, data));
    }

    /** \brief Number of data sets in the container
     **/
    size_t size() const { return m_sets.size(); }

    /** \brief STL-style iterator to first data set
     **/
    iterator begin() const {
        return m_sets.begin();
    }

    /** \brief STL-style iterator to end
     **/
    iterator end() const {
        return m_sets.end();
    }

    /** \brief Access to the name of the current data set (const)
     **/
    const std::string &name(iterator i) const {
        return i->name;
    }

    /** \brief Access to the current data set (const)
     **/
    const arma::vec &data(iterator i) const {
        return i->data;
    }

    /** \brief Print population data as table
        \param out Output stream
        \param l Row labels
        \param colwidth Column width (default: 16)
        \param prec Precision (default: 6)
        \param offset Line offset (default: 2)
     **/
    void print(std::ostream &out, const std::vector<std::string> &l,
        size_t colwidth = 16, size_t prec = 6, size_t offset = 2) const;
};

} // namespace libwfa

#endif // LIBWFA_POP_DATA_H
