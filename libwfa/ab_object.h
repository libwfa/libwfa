#ifndef LIBWFA_AB_OBJECT_H
#define LIBWFA_AB_OBJECT_H


namespace libwfa {

/** \brief Container for objects with spin property
    \tparam T Object type

    The container stores two objects with alpha- and beta-spin and a flag to
    indicate, if both objects have to be identical. The case of alpha == beta
    will be enforced by the container by using the same object for alpha-spin
    and beta-spin.

    The object type has to possess default and copy constructor for the class
    to work properly.

    \ingroup libwfa
 **/
template<typename T>
class ab_object {
private:
    bool m_aeqb; //!< If alpha-spin equals beta spin
    T* m_data_a;
    T* m_data_b;

public:
    ab_object(bool aeqb = false) : m_aeqb(aeqb) {
        m_data_a = new T();
        m_data_b = (m_aeqb ? m_data_a : new T());
    }

    ~ab_object() {
        delete m_data_a;
        if (! m_aeqb) delete m_data_b;

        m_data_a = m_data_b = 0;
    }

    /** \brief Set alpha == beta

        Modifies the container to enforce both objects to be identical. The
        contents of the beta objects is deleted and only the alpha object is
        kept.
     **/
    void set_alpha_eq_beta() {
        if (m_aeqb) return;

        delete m_data_b;
        m_data_b = m_data_a;
        m_aeqb = true;
    }

    /** \brief Set alpha != beta

        Modifies the container to allow the objects with alpha and beta-spin
        to change independently. However, both objects will be identical
        copies of one another first.
     **/
    void set_alpha_neq_beta() {
        if (! m_aeqb) return;

        m_data_b = new T(*m_data_a);
        m_aeqb = false;
    }

    /** \brief Are alpha- and beta-spin matrices identical
     **/
    bool is_alpha_eq_beta() const { return m_aeqb; }

    /** \brief Return alpha-spin matrix
     **/
    T &alpha() { return *m_data_a; }

    /** \brief Return alpha-spin matrix (const version)
     **/
    const T &alpha() const { return *m_data_a; }

    /** \brief Return beta-spin matrix
     **/
    T &beta() { return *m_data_b; }

    /** \brief Return beta-spin matrix (const version
     **/
    const T &beta() const { return *m_data_b; }
};

} // namespace adcman

#endif // LIBWFA_AB_MATRIX_H

