#ifndef LIBWFA_CTNUM_DATA_H
#define LIBWFA_CTNUM_DATA_H

#include <armadillo>
#include <list>
#include <string>

namespace libwfa {

/** \brief Container for several named population data sets

    \ingroup libwfa
 **/
class ctnum_data {
private:
    /** \brief Structure to keep one single data set
     **/
    struct data_set{
        std::string name;
        double energy;
        double osc_strength;
        arma::Mat<double> data;

        data_set(const std::string &name_, double ene, double osc,
            const arma::Mat<double> &data_ = arma::Mat<double>()) :
            name(name_), energy(ene), osc_strength(osc), data(data_) {
        }
    };

public:
    typedef std::list<data_set>::const_iterator iterator; //!< Iterator to traverse the data

private:
    std::list<data_set> m_sets; //!< Collection of data sets

public:
    /** \brief Add a new empty set for CT number data
        \param name Name of data set
        \param ene Excitation energy
        \param osc Oscillator strength
        \return Reference to the empty data set.
     **/
    arma::Mat<double> &add(const std::string &name, double ene, double osc) {
        m_sets.push_back(data_set(name, ene, osc));
        return m_sets.back().data;
    }

    /** \brief Add a new set for CT number data
        \param name Name of data set
        \param ene Excitation energy
        \param osc Oscillator strength
        \param data The data to be added (copied)
     **/
    void add(const std::string &name, double ene, double osc,
            const arma::Mat<double> &data) {
        m_sets.push_back(data_set(name, ene, osc, data));
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

    /** \brief Access to the excitation energy associated with the current
            data set (const)
     **/
    double energy(iterator i) const {
        return i->energy;
    }

    /** \brief Access to the oscillator strength associated with the current
            data set (const)
     **/
    double osc_strength(iterator i) const {
        return i->osc_strength;
    }

    /** \brief Access to the current data set (const)
     **/
    const arma::Mat<double> &data(iterator i) const {
        return i->data;
    }
};

} // namespace libwfa

#endif // LIBWFA_CTNUM_DATA_H
