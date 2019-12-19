#include <iomanip>
#include <libwfa/libwfa_exception.h>
#include "cube_writer.h"
#include "export_cube_base.h"

using namespace arma;

namespace libwfa {


const char export_cube_base::k_clazz[] = "export_cube_base";


export_cube_base::export_cube_base(const grid3d &grid,
    const uvec &atnum, const mat &coords,
    const std::string &prefix, size_t batchsz, size_t nmax,
    const std::string &comment) :
    m_grid(grid), m_atnum(atnum), m_coords(coords), m_prefix(prefix),
    m_comment(comment), m_batchsz(batchsz), m_nmax(nmax) {

    if (m_batchsz == 0) m_batchsz = m_grid.size();
    if (m_grid.size() < m_batchsz) m_batchsz = m_grid.size();
}


void export_cube_base::perform(const std::string &name,
    const std::string &desc, const mat &data, int do_esp) {

    static const char method[] = "add(const std::string &, "
            "const std::string &, const mat &)";

    if (data.n_rows != data.n_cols) {
        libwfa_exception(k_clazz, method, __FILE__, __LINE__, "data");
    }
    if (m_dms.find(name) != m_dms.end()) {
        libwfa_exception(k_clazz, method, __FILE__, __LINE__, name.c_str());
    }

    dm_data *ptr = new dm_data(desc, data);
    m_dms.insert(dm_list::value_type(name, ptr));

    if (m_nmax != 0 && m_dms.size() + m_orbs.size() > m_nmax) do_export();

    if (do_esp>0) {
        export_esp(name, desc, data, do_esp);
    }
}


void export_cube_base::perform(const std::string &name,
    const std::string &desc, const std::vector<size_t> &idx,
    const mat &vecs) {

    static const char method[] =
            "add(const std::string &, const std::string &, "
            "const std::vector<size_t> &, const mat &)";

    if (idx.size() != vecs.n_cols) {
        libwfa_exception(k_clazz, method, __FILE__, __LINE__, "idx");
    }
    if (m_orbs.find(name) != m_orbs.end()) {
        libwfa_exception(k_clazz, method, __FILE__, __LINE__, name.c_str());
    }

    orb_data *ptr = new orb_data(desc, idx, vecs);
    m_orbs.insert(orb_list::value_type(name, ptr));

    if (m_nmax != 0 && m_dms.size() + m_orbs.size() > m_nmax) do_export();
}


void export_cube_base::do_export() {

    typedef std::vector<cube_writer *> cube_writer_list;

    if (m_dms.size() + m_orbs.size() == 0) return;

    // Prepare list of output files for dms
    cube_writer_list dm_writers;
    for (dm_list::iterator i = m_dms.begin(); i != m_dms.end(); i++) {
        std::string fname(m_prefix + i->first + ".cube");
        dm_writers.push_back(new cube_writer(fname, i->second->desc,
                m_comment, m_grid, m_atnum, m_coords));
    }

    // Prepare list of output files for orbitals
    size_t ndigits = 0, nao = 0;
    if (m_orbs.size() != 0) nao = m_orbs.begin()->second->data.n_rows;
    else nao = m_dms.begin()->second->data.n_rows;
    for (size_t nn = nao; nn > 0; nn /= 10, ndigits++) { ; }

    cube_writer_list orb_writers;
    for (orb_list::iterator i = m_orbs.begin(); i != m_orbs.end(); i++) {

        orb_data &data = *(i->second);
        for (std::vector<size_t>::iterator j = data.idx.begin();
                j != data.idx.end(); j++) {

            std::ostringstream ss1, ss2;
            ss1 << m_prefix << i->first << "." << std::setw(ndigits) <<
                    std::setfill('0') << *j << ".cube";
            ss2 << data.desc << " (orbital " << *j << ")";
            orb_writers.push_back(new cube_writer(ss1.str(), ss2.str(),
                    m_comment, m_grid, m_atnum, m_coords));
        }
    }

    // Now perform the loop over the batches
    mat pts(3, m_batchsz, fill::zeros);
    mat b2g(m_batchsz, nao, fill::zeros);

    size_t ipts = 0;
    while (ipts < m_grid.size()) {

        // 1) Build grid points
        size_t npts = m_grid.build_pts(ipts, pts);

        // 2) Precompute the basis functions on grid points
        //    FP: This is independent from the density matrices
        evaluate_on_grid(pts, npts, b2g);

        // 3) Loop over density matrix data
        //    FP: This is a simple matrix multiplication
        cube_writer_list::iterator id = dm_writers.begin();
        for (dm_list::iterator i = m_dms.begin();
                i != m_dms.end(); i++, id++) {

            vec data(npts);
            data = sum((b2g * i->second->data) % b2g, 1);
            (*id)->write(data);
        }

        // 4) Loop over orbital data
        cube_writer_list::iterator io = orb_writers.begin();
        for (orb_list::iterator i = m_orbs.begin();
                i != m_orbs.end(); i++) {

            mat data = b2g.rows(0, npts - 1) * i->second->data;
            for (size_t j = 0; j < data.n_cols; j++, io++)
                (*io)->write(data.col(j));
        }

        ipts += npts;
    }

    // Delete cube writers and close files
    for (cube_writer_list::iterator i = dm_writers.begin();
            i != dm_writers.end(); i++) {
        delete *i; *i = 0;
    }
    dm_writers.clear();

    for (cube_writer_list::iterator i = orb_writers.begin();
            i != orb_writers.end(); i++) {
        delete *i; *i = 0;
    }
    orb_writers.clear();

    // Delete density matrix data and orbital data
    clear_data();
}

void export_cube_base::export_esp(const std::string &name, const std::string &desc,
        const arma::mat &dens, int do_esp) {

    mat pts(3, m_batchsz, fill::zeros);

    std::string fname(m_prefix + name + "_esp.cube");
    cube_writer espw(fname, desc, m_comment, m_grid, m_atnum, m_coords);

    size_t ipts = 0;
    while (ipts < m_grid.size()) {
        size_t npts = m_grid.build_pts(ipts, pts);

        vec esp(npts);
        evaluate_esp(pts, npts, dens, esp);
        if (do_esp==2) {
            nuc_esp(pts, npts, esp);
        }
        espw.write(esp);

        ipts += npts;
    }
}

void export_cube_base::nuc_esp(const arma::mat &pts, size_t npts, arma::vec &esp) {
    // The nuclear contribution to the ESP.
    // This could be stored somewhere. But it probably does not hurt to recompute
    //  it for every density.

    // Loop over grid points
    for (size_t k = 0; k < npts; k++) {
        // Loop over atoms
        for (size_t iat = 0; iat < m_atnum.size(); iat++) {
            esp(k) += m_atnum(iat) / arma::norm(pts.col(k) - m_coords.col(iat));
        }
    }
}

void export_cube_base::clear_data() {

    for (dm_list::iterator i = m_dms.begin(); i != m_dms.end(); i++) {
        delete i->second; i->second = 0;
    }
    m_dms.clear();

    for (orb_list::iterator i = m_orbs.begin(); i != m_orbs.end(); i++) {
        delete i->second; i->second = 0;
    }
    m_orbs.clear();
}


} // namespace libwfa
