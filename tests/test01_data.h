#ifndef LIBWFA_TEST01_DATA_H
#define LIBWFA_TEST01_DATA_H

#include "test_data_base.h"

namespace libwfa {


/** \brief Data for test system 1

    The test system is He-Li at a distance of 4 Angstrom computed using an
    STO-3G basis set and ADC(2)-s as excited state method.

    The following input file was used to generate the test data
    \code
    $rem
    jobtype = sp
    method = adc(2)
    basis = sto-3g
    ee_singlets = 4
    cc_symmetry = false
    adc_prop_es = true
    adc_nguess_singles = 4
    adc_nguess_doubles = 4
    adc_davidson_conv = 6
    make_cube_files true
    $end

    $molecule
    0 2
    He
    Li 1  4.000
    $end

    $plots
    Plot all four types of quantities
    80 -6.0 6.0
    40 -3.0 3.0
    40 -3.0 3.0
    0 0 2 2
    1 2
    1 2
    $end
    \endcode

    \ingroup libwfa_tests
 **/
class test01_data : public test_data_base {
public:
    static const size_t k_nao; //!< Number of AOs
    static const size_t k_nmo; //!< Number of MOs
    static const size_t k_natoms; //!< Number of atoms
    static const size_t k_nstat; //!< Number of states    

private:
    arma::uvec m_atnum; //!< Atomic numbers
    arma::vec m_nch; //!< Nuclear charges
    arma::uvec m_bf2nuc; //!< Map of basis functions to nuclei

    arma::mat m_pop_mulliken[2]; //!< Reference for population analysis
    arma::mat m_pop_loewdin[2]; //!< Reference for population analysis

public:
    test01_data();

    bool aeqb() { return false; }

    size_t nstates() { return k_nstat; }

    arma::uvec atomic_numbers() { return m_atnum; }

    arma::vec nuclear_charges() { return m_nch; }

    arma::uvec bf2nuclei() { return m_bf2nuc; }

    arma::vec pop_mulliken(size_t istate, bool alpha) {
        if (alpha) return m_pop_mulliken[0].col(istate);
        else return m_pop_mulliken[1].col(istate);
    }

    arma::vec pop_loewdin(size_t istate, bool alpha) {
        if (alpha) return m_pop_loewdin[0].col(istate);
        else return m_pop_loewdin[1].col(istate);
    }
};


} // namespace libwfa

#endif // LIBWFA_TEST01_DATA_H
