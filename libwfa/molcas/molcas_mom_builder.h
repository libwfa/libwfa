#ifndef LIBWFA_MOLCAS_MOM_BUILDER_H
#define LIBWFA_MOLCAS_MOM_BUILDER_H

#include "H5Cpp.h"
#include <libwfa/core/mom_builder_i.h>

namespace libwfa {

/** \brief Implementation of mom_builder_i for Molcas
  
   Store arrays in original symmetrized form?

    \ingroup libwfa
 **/
class molcas_mom_builder : public libwfa::mom_builder_i {
private:
    mutable std::vector<arma::mat *> m_op;
    unsigned int m_maxmm;
    size_t m_nao;
    /*H5::File m_file; //!< HDF5 file with data
    const arma::mat  m_desym; //!< Desymmetrization matrix*/

public:
    /**\brief Constructor
     **/
    molcas_mom_builder(size_t nao, unsigned int maxmm) :
        m_nao(nao), m_maxmm(maxmm), m_op(maxmm * 3 + 1, 0) { }

    virtual ~molcas_mom_builder() {}

    /** \brief Set the desired operator matrix
        \param c Coordinate direction
        \param n Exponent
        \return Operator matrix
     **/
    arma::mat &set(size_t c, size_t n);
    
    /** \brief Returns the desired operator as a matrix
        \param c Coordinate direction
        \param n Exponent
        \return Operator matrix
     **/
    const arma::mat &get(size_t c, size_t n) const;

    /** \copydoc mom_builder_i::max_moment
     **/
    virtual size_t max_moment() const { return m_maxmm; }

    /** \copydoc mom_builder_i::perform
     **/
    virtual double perform(const arma::mat &dm,
        size_t c1, size_t n1, size_t c2, size_t n2) const;

    /** \copydoc mom_builder_i::perform
     **/
    virtual double perform(const arma::mat &dm, size_t c, size_t n) const;

private:
    /** \brief Returns the index
        \param c Coordinate direction
        \param n Exponent
        \return Index
     **/    
    static size_t get_index(size_t c, size_t n) {
        if (n == 0) return 0;
        else return n * 3 - 2 + c;
    }

    void initialize(size_t pos) const; // -> delete
};


} // namespace libwfa

#endif // LIBWFA_MOLCAS_MOM_BUILDER_H
