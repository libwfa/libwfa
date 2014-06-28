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

public:
    test01_data() : test_data_base("test01") { }
};


} // namespace libwfa

#endif // LIBWFA_TRANSFORMATIONS_DM_TEST_H
