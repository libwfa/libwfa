#ifndef LIBWFA_TEST01_DATA_H
#define LIBWFA_TEST01_DATA_H

#include "test_data_base.h"

namespace libwfa {


/** \brief Data for test system 1

    The test system is He-Li at a distance of 4 Angstrom computed using an
    STO-3G basis set and ADC(2)-s as excited state method.

    \ingroup libwfa_tests
 **/
class test01_data : public test_data_base {
public:
    static const size_t k_nao; //!< Number of AOs
    static const size_t k_nmo; //!< Number of MOs
    static const size_t k_natoms; //!< Number of atoms

private:
    arma::Col<size_t> m_atnum; //!< Atomic numbers
    arma::Col<double> m_nch; //!< Nuclear charges
    arma::Col<size_t> m_bf2nuc; //!< Map of basis functions to nuclei

public:
    test01_data();

    bool aeqb() { return false; }

    size_t nstates() { return 2; }

    arma::Col<size_t> atomic_numbers() { return m_atnum; }

    arma::Col<double> nuclear_charges() { return m_nch; }

    arma::Col<size_t> bf2nuclei() { return m_bf2nuc; }
};


} // namespace libwfa

#endif // LIBWFA_TEST01_DATA_H
