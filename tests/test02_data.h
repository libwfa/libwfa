#ifndef LIBWFA_TEST02_DATA_H
#define LIBWFA_TEST02_DATA_H

#include "test_data_base.h"

namespace libwfa {


/** \brief Data for test system 2

    The test system is H2O computed using 3-21g basis set and ADC(2)-s as
    excited state method.

    \ingroup libwfa_tests
 **/
class test02_data : public test_data_base {
public:
    static const size_t k_nao; //!< Number of AOs
    static const size_t k_nmo; //!< Number of MOs
    static const size_t k_natoms; //!< Number of atoms

private:
    arma::Col<size_t> m_atnum; //!< Atomic numbers
    arma::Col<double> m_nch; //!< Nuclear charges
    arma::Col<size_t> m_bf2nuc; //!< Map of basis functions to nuclei

public:
    test02_data();

    bool aeqb() { return true; }

    size_t nstates() { return 2; }

    arma::Col<size_t> atomic_numbers() { return m_atnum; }

    arma::Col<double> nuclear_charges() { return m_nch; }

    arma::Col<size_t> bf2nuclei() { return m_bf2nuc; }
};


} // namespace libwfa

#endif // LIBWFA_TEST02_DATA_H
