#ifndef EX_ANA_PRI_TEST_H_
#define EX_ANA_PRI_TEST_H_
#include <libtest/unit_test.h>

namespace libwfa {

/** \brief Tests for libwfa::ex_analyse class

    \ingroup libwfa_tests
 **/
class ex_ana_pri_test: public libtest::unit_test {
public:
    /** \brief Tests printer and saves output in ex_ana_pri_test.txt
     **/
    virtual void perform() throw(libtest::test_exception);
private:

}; //end class


}



#endif /* EX_ANA_PRI_TEST_H_ */
