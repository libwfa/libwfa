#ifndef LIBWFA_AB_OBJECT_H
#define LIBWFA_AB_OBJECT_H


namespace libwfa {

/** \brief Container for objects with spin property
    \tparam T Object type

    The container stores two objects with alpha- and beta-spin and a flag to
    indicate, if both objects have to be identical. The case of alpha == beta
    will be enforced by the container by using the same object for alpha-spin
    and beta-spin.

    The object type has to possess default and copy constructor to work properly.

    \ingroup libwfa
 **/
template<typename T>
class ab_object {
private:
    bool m_aeqb; //!< If alpha-spin equals beta spin
    T* m_data_a; //!< Object data for alpha spin
    T* m_data_b; //!< Object data for beta spin

public:
    /** \brief Default constructor
        \param aeqb Is alpha equal to beta (default: false)
     **/
    ab_object(bool aeqb = false) : m_aeqb(aeqb) {
        m_data_a = new T();
        m_data_b = (m_aeqb ? m_data_a : new T());
    }

    /** \brief Copy constructor
        \param other Object to copy
        
        Performs deep copy
     **/
    ab_object(const ab_object &other);

    /** \brief Destructor
     **/
    ~ab_object() {
        delete m_data_a;
        if (! m_aeqb) delete m_data_b;

        m_data_a = m_data_b = 0;
    }

    /** \brief Assignment operator
        \param other Object to get data from
     **/
    ab_object &operator=(const ab_object &other);
    
    /** \brief Assignment of ab_object
        \param other Object to get data from
     **/
    void set(const ab_object &other) {
        if (m_aeqb && ! other.is_alpha_eq_beta()) {
            m_data_b = new T(other.beta());
        }
        else if (! m_aeqb && other.is_alpha_eq_beta()) {
            delete m_data_b;
            m_data_b = m_data_a;
        }
        else if (! m_aeqb && ! other.is_alpha_eq_beta()) {
            *m_data_b = other.beta();
        }
        m_aeqb = other.is_alpha_eq_beta();
        *m_data_a = other.alpha();
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


template<typename T>
ab_object<T>::ab_object(const ab_object<T> &other) :
    m_aeqb(other.is_alpha_eq_beta()) {

    m_data_a = new T(other.alpha());
    m_data_b = (m_aeqb ? m_data_a : new T(other.beta()));
}


template<typename T>
ab_object<T> &ab_object<T>::operator=(const ab_object<T> &other) {

    if (m_aeqb) {
        if (! other.is_alpha_eq_beta()) m_data_b = new T(other.beta());
    }
    else {
        if (other.is_alpha_eq_beta()) {
            delete m_data_b; m_data_b = m_data_a;
        }
        else { *m_data_b = other.beta(); }
    }
    *m_data_a = other.alpha();
    m_aeqb = other.is_alpha_eq_beta();
}


} // namespace libwfa

#endif // LIBWFA_AB_MATRIX_H

