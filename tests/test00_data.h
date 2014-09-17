#ifndef LIBWFA_TEST00_DATA_H
#define LIBWFA_TEST00_DATA_H

#include "test_data_base.h"

namespace libwfa {


/** \brief Random test data

    This test system comprises only overlap matrix and MO coefficients
    constructed from random data.

    \ingroup libwfa_tests
 **/
class test00_data : public test_data_base {
public:
    static const size_t k_nao; //!< Number of AOs
    static const size_t k_nmo; //!< Number of MOs

public:
    test00_data();

    bool aeqb() { return false; }

};


} // namespace libwfa

#endif // LIBWFA_TEST00_DATA_H
