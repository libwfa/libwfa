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
        arma::Col<double> data;

        data_set(const std::string &name_,
            const arma::Col<double> &data_ = arma::Col<double>()) :
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
    arma::Col<double> &add(const std::string &name) {
        m_sets.push_back(data_set(name));
        return m_sets.back().data;
    }

    /** \brief Add a new set for population data
        \param name Name of data set
        \param data The data to be added
     **/
    void add(const std::string &name, const arma::Col<double> &data) {
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
    const arma::Col<double> &data(iterator i) const {
        return i->data;
    }
};

} // namespace libwfa

#endif // LIBWFA_POP_DATA_H
