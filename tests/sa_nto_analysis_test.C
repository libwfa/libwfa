#include <sstream>
#include <libwfa/analyses/sa_nto_analysis.h>
#include <libwfa/libwfa_exception.h>
#include "sa_nto_analysis_test.h"
#include "test01_data.h"
#include "test02_data.h"

namespace libwfa {

using namespace arma;


void sa_nto_analysis_test::perform() throw(libtest::test_exception) {

    test_1<test01_data>();
    test_1<test02_data>();
    //fail_test("sa_nto_analysis_test::perform()", __FILE__, __LINE__, "NIY");
}

template<typename TestData>
void sa_nto_analysis_test::test_1() throw(libtest::test_exception) {

    static const char *testname = "sa_nto_analysis_test::test_1()";

    try {

    // Prepare input data:

    size_t nao = TestData::k_nao;
    size_t nmo = TestData::k_nmo;
    TestData data;

    mat s(nao, nao);
    read_matrix(data, testname, "s", s);

    ab_matrix c(data.aeqb());
    c.alpha() = mat(nao, nmo);
    if (! data.aeqb()) c.beta() = mat(nao, nmo);
    read_ab_matrix(data, testname, "c", c);

    ab_matrix hdm_av(data.aeqb()), edm_av(data.aeqb());
    hdm_av.alpha() = mat(nao, nao, fill::zeros);
    edm_av.alpha() = mat(nao, nao, fill::zeros);
    if (! data.aeqb()) {
        hdm_av.beta() = mat(nao, nao, fill::zeros);
        edm_av.beta() = mat(nao, nao, fill::zeros);
    }

    // Form state-averaged hole and particle densities
    for (size_t istate = 1; istate <= data.nstates(); istate++) {
        ab_matrix tdm(data.aeqb());
        tdm.alpha() = mat(nao, nao);
        if (! data.aeqb()) tdm.beta() = mat(nao, nao);
        std::ostringstream ssdm; ssdm << "tdm" << istate;
        read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);

        ab_matrix hdm, edm;
        nto_analysis::form_eh(s, tdm.t(), edm, hdm);

        hdm_av += hdm;
        edm_av += edm;
    }

    // Form spin-traced hole and particle densities
    ab_matrix edm_av_s(true), hdm_av_s(true);
    edm_av_s.alpha() = 0.5 * (edm_av.alpha() + edm_av.beta());
    hdm_av_s.alpha() = 0.5 * (hdm_av.alpha() + hdm_av.beta());

    // Construct SA-NTOs

    sa_nto_analysis sa_nto1(s, nto_analysis(s, c, edm_av, hdm_av));
    sa_nto_analysis sa_nto2(s, nto_analysis(s, c, edm_av_s, hdm_av_s));

    const ab_matrix &u1 = sa_nto1.get_transf_l();
    const ab_matrix &v1 = sa_nto1.get_transf_r();
    const ab_matrix &u2 = sa_nto2.get_transf_l();
    const ab_matrix &v2 = sa_nto2.get_transf_r();

    if (accu(abs(u1.alpha().t() * u1.alpha() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "U(1,alpha) not unitary.");
    if (accu(abs(u1.beta().t() * u1.beta() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "U(1,beta) not unitary.");
    if (accu(abs(v1.alpha() * v1.alpha().t() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "V(1,alpha) not unitary.");
    if (accu(abs(v1.beta() * v1.beta().t() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "V(1,beta) not unitary.");
    if (! u2.is_alpha_eq_beta() || ! v2.is_alpha_eq_beta())
        fail_test(testname, __FILE__, __LINE__, "SA-NTO(2): Alpha != beta.");
    if (accu(abs(u2.alpha().t() * u2.alpha() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "U(2,alpha) not unitary.");
    if (accu(abs(u2.beta().t() * u2.beta() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "U(2,beta) not unitary.");
    if (accu(abs(v2.alpha() * v2.alpha().t() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "V(2,alpha) not unitary.");
    if (accu(abs(v2.beta() * v2.beta().t() - s) > 1e-12) != 0)
        fail_test(testname, __FILE__, __LINE__, "V(1,beta) not unitary.");

    // Loop to decompose TDM for each state

    for (size_t istate = 1; istate <= data.nstates(); istate++) {

        ab_matrix tdm(data.aeqb());
        tdm.alpha() = mat(nao, nao);
        if (! data.aeqb()) tdm.beta() = mat(nao, nao);
        std::ostringstream ssdm; ssdm << "tdm" << istate;
        read_ab_matrix(data, testname, ssdm.str().c_str(), tdm);

        tdm = tdm.t();

        sa_nto_analysis sa_nto0(s, nto_analysis(s, c, tdm));

        ab_matrix x0, x1, x2;
        sa_nto0.decompose(tdm, x0);
        sa_nto1.decompose(tdm, x1);
        sa_nto2.decompose(tdm, x2);

        const ab_matrix &ui = sa_nto0.get_transf_l();
        const ab_matrix &vi = sa_nto0.get_transf_r();
        if (accu(abs(ui.alpha().t() * ui.alpha() - s) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "U(0, alpha) not unitary.");
        if (accu(abs(ui.beta().t() * ui.beta() - s) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "U(0, beta) not unitary.");
        if (accu(abs(vi.alpha() * vi.alpha().t() - s) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "V(0, alpha) not unitary.");
        if (accu(abs(vi.beta() * vi.beta().t() - s) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "V(0, beta) not unitary.");
        if (accu(abs(ui.alpha() * tdm.alpha() * vi.alpha() - x0.alpha()) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform (alpha).");
        if (accu(abs(ui.beta() * tdm.beta() * vi.beta() - x0.beta()) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform (beta).");
        if (accu(abs(x0.alpha() - diagmat(x0.alpha().diag()))) > 1e-12)
            fail_test(testname, __FILE__, __LINE__, "Not diagonal (alpha)");
        if (accu(abs(x0.beta() - diagmat(x0.beta().diag()))) > 1e-12)
            fail_test(testname, __FILE__, __LINE__, "Not diagonal (beta)");
        if (accu(abs(u1.alpha() * tdm.alpha() * v1.alpha() - x1.alpha()) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform (1, alpha).");
        if (accu(abs(u1.beta() * tdm.beta() * v1.beta() - x1.beta()) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform (1, beta).");
        if (accu(abs(u2.alpha() * tdm.alpha() * v2.alpha() - x2.alpha()) > 1e-12) != 0) {
            fail_test(testname, __FILE__, __LINE__, "Bad transform (2, alpha).");
        }
        if (accu(abs(u2.beta() * tdm.beta() * v2.beta() - x2.beta()) > 1e-12) != 0)
            fail_test(testname, __FILE__, __LINE__, "Bad transform (2, beta).");
    }

    } catch(libtest::test_exception &e) {
        throw;
    } catch(std::exception &e) {
        fail_test(testname, __FILE__, __LINE__, e.what());
    }

}


} // namespace libwfa
