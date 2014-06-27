#ifndef LIBWFA_TEST_DATA_BASE_H
#define LIBWFA_TEST_DATA_BASE_H

#include <armadillo>
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
    bool read_matrix(const char *testname,
        const char *fname, arma::Mat<double> &m);

    /** \brief Read double from file
        \param[in] testname Name of the test
        \param[in] fname File name to read from (w/o prefix)
        \param[out] d Variable to put data
        \return True if successful
     **/
    bool read_double(const char *testname, const char *fname, double &d);

private:
    std::string make_filename(const char *fname) const;
};


} // namespace libwfa

#endif // LIBWFA_TRANSFORMATIONS_DM_TEST_H
