#ifndef LIBWFA_TEST_DATA_BASE_H
#define LIBWFA_TEST_DATA_BASE_H

#include <armadillo>
#include <libwfa/core/ab_matrix.h>
#include <string>

namespace libwfa {


/** \brief Base class for test data

    \ingroup libwfa_tests
 **/
class test_data_base {
private:
    std::string m_prefix; //!< Prefix to append to each filename

public:
    /** \brief Constructor
        \param prefix File prefix for test
     **/
    test_data_base(const std::string &prefix) : m_prefix(prefix) { }

    /** \brief Virtual destructor
     **/
    virtual ~test_data_base() { }

    /** \brief Read armadillo matrix or vector from file
        \param[in] testname Name of the test
        \param[in] fname File name to read from (w/o prefix)
        \param[out] m Matrix to put data read
        \return True if successful
     **/
    template<typename T>
    bool read_matrix(const char *testname,
        const char *fname, arma::Mat<T> &m) const;

    /** \brief Read double from file
        \param[in] testname Name of the test
        \param[in] fname File name to read from (w/o prefix)
        \param[out] d Variable to put data
        \return True if successful
     **/
    bool read_double(const char *testname,
        const char *fname, double &d) const;

    /** \brief Construct the complete filename from prefix and basic filename
        \param fname Basic filename
        \return Complete filenam
     **/
    std::string make_filename(const char *fname) const;
};


template<typename T>
bool test_data_base::read_matrix(const char *testname,
    const char *fname, arma::Mat<T> &m) const {

    std::ifstream in(make_filename(fname).c_str());
    if (in.fail()) {
        std::cout << "in.fail()" << std::endl;
        return false;
    }

    int nr = 0, nc = 0;
    in >> nr  >> nc;

    if (nr != m.n_rows) {
        std::cout << "Wrong number of rows" << std::endl;
        return false;
    }
    if (nc != m.n_cols) {
        std::cout << "Wrong number of cols" << std::endl;
        return false;
    }

    T *ptr = m.memptr();

    for (size_t i = 0; i < m.n_elem && in.good(); i++, ptr++) in >> *ptr;

    if (in.fail()) {
        std::cout << "in.fail() 2" << std::endl;
        return false;
    }

    return true;
}


} // namespace libwfa

#endif // LIBWFA_TEST_DATA_BASE_H
