#ifndef LIBWFA_MOLDEN_FILE_BASE_H
#define LIBWFA_MOLDEN_FILE_BASE_H

#include "ab_matrix.h"
#include "ab_vector.h"

namespace libwfa {

namespace cube {

/** \brief Grid in real space (3D)

    \ingroup libwfa
 **/
struct grid3d {
    double ranges[6]; //!< [min, max] ranges in each direction
    unsigned int npts[3]; //!< Number of grid points in each direction

    /** \brief Constructor
     **/
    grid3d();

    //! \name Setters (with error checks)
    //@{
    void set_xrange(double min, double max, unsigned int npts);
    void set_yrange(double min, double max, unsigned int npts);
    void set_zrange(double min, double max, unsigned int npts);
    //@}

    /** \brief Consistency check
     **/
    void check() const;
};


/** \brief Data type to generate cube data

    \ingroup libwfa
 **/
class data_type {
private:
    char m_type; //!< Data type

private:
    explicit data_type(char type) : m_type(type) { }

public:
    bool operator==(const data_type &t) const { return m_type == t.m_type; }
    bool operator!=(const data_type &t) const { return m_type != t.m_type; }

public:
    static const data_type orb_a; //!< Alpha orbital coefficients
    static const data_type orb_b; //!< Beta orbital coefficients
    static const data_type sdm_a; //!< Alpha state density matrix
    static const data_type sdm_b; //!< Beta state density matrix
    static const data_type tdm_a; //!< Alpha transition density matrix
    static const data_type tdm_b; //!< Beta transition density matrix
    static const data_type adm_a; //!< Alpha transition density matrix
    static const data_type adm_b; //!< Beta transition density matrix
    static const data_type ddm_a; //!< Alpha transition density matrix
    static const data_type ddm_b; //!< Beta transition density matrix
};

} // namespace cube

/** \brief Base class to export data as cube files.

    Base class to export certain data (see data_type for possible data) as
    grid data into cube files.

    \ingroup libwfa
 **/
class export_cube_base {
private:
    cube::grid3d m_grid;

public:
    /** \brief Constructor
        \param g Grid of the cube
     **/
    export_cube_base(const cube::grid3d &g) : m_grid(g) { }

    virtual ~export_cube_base() { }

    /** \brief Return the grid information
     **/
    const cube::grid3d &grid() const { return m_grid; }

    /** \brief Export several sets of data provided as matrix
        \param type Type of data
        \param idx Indices for each data
        \param data Matrix containing the data objects
     **/
    virtual void perform(cube::data_type type, const std::vector<size_t> &idx,
            const arma::Mat<double> &data) = 0;
};

} // namespace adcman

#endif
