#ifndef LIBWFA_TEST01_DATA_H
#define LIBWFA_TEST01_DATA_H

#include "test_data_base.h"

namespace libwfa {


/** \brief

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
