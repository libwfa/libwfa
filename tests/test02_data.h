#ifndef LIBWFA_TEST02_DATA_H
#define LIBWFA_TEST02_DATA_H

#include "test_data_base.h"

namespace libwfa {


/** \brief Data for test system 2

    The test system is H2O computed using 3-21g basis set and ADC(2)-s as
    excited state method. No point-group symmetry has been used. Two states
    have been computed (the lowest 1 singlet and the lowest triplet state).

    The following input file was used to generate the test data
    \code
    $rem
    jobtype = sp
    method = adc(2)
    basis = 3-21g
    ee_singlets = 1
    ee_triplets = 1
    n_frozen_core = fc
    cc_symmetry = false
    adc_prop_es = true
    adc_nguess_singles = 2
    adc_nguess_doubles = 1
    adc_davidson_conv = 6
    adc_davidson_thresh = 8
    make_cube_files true
    $end

    $molecule
    0 1
    O
    H  1  0.957
    H  1  0.957  2   104.5
    $end

    $plots
    Plot all four types of quantities
    40 -3.0 3.0
    40 -3.0 3.0
    40 -3.0 3.0
    0 0 2 2
    1 2
    1 2
    $end
    \endcode

    \ingroup libwfa_tests
 **/
class test02_data : public test_data_base {
public:
    static const size_t k_nao; //!< Number of AOs
    static const size_t k_nmo; //!< Number of MOs
    static const size_t k_natoms; //!< Number of atoms
    static const size_t k_nstat; //!< Number of states    

private:
    arma::Col<size_t> m_atnum; //!< Atomic numbers
    arma::Col<double> m_nch; //!< Nuclear charges
    arma::Col<size_t> m_bf2nuc; //!< Map of basis functions to nuclei
    
    arma::Mat<double> m_popref_a; //<! Reference for population analysis

public:
    test02_data();

    bool aeqb() { return true; }

    size_t nstates() { return k_nstat; }

    arma::Col<size_t> atomic_numbers() { return m_atnum; }

    arma::Col<double> nuclear_charges() { return m_nch; }

    arma::Col<size_t> bf2nuclei() { return m_bf2nuc; }
    
    arma::Col<double> popref(size_t istate, bool alpha) {
        return m_popref_a.col(istate);
        // this test has no beta field
    }
};


} // namespace libwfa

#endif // LIBWFA_TEST02_DATA_H
