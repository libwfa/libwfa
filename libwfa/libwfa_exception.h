#ifndef LIBWFA_LIBWFA_EXCEPTION_H
#define LIBWFA_LIBWFA_EXCEPTION_H

#include <exception>
#include <memory>

namespace libwfa {


/**	\brief Base exception class

	\ingroup libutil_exceptions
 **/
class libwfa_exception : public std::exception {
private:
    char m_clazz[128]; //!< Class name
    char m_method[128]; //!< Method name
    char m_file[128]; //!< Source file name
    unsigned int m_line; //!< Line number
    char m_message[256]; //!< Exception message
    char m_what[1024]; //!< Composed message available via what()

public:
	//!	\name Construction and destruction
	//@{

    /** \brief Creates an %exception using full details
        \param clazz Class name.
        \param method Method name.
        \param file Source file name.
        \param line Line number in the source file.
        \param message Error message.
     **/
    libwfa_exception(const char *clazz, const char *method,
        const char *file, unsigned int line, const char *message) throw() {

        init(clazz, method, file, line, message);
    }

    /** \brief Copy constructor
        \param e Other exception
     **/
    libwfa_exception(const libwfa_exception &e) {
        init(e.m_clazz, e.m_method, e.m_file, e.m_line, e.m_message);
    }

    /**	\brief Virtual destructor
	 **/
	virtual ~libwfa_exception() throw() { };

	//@}


	//!	\name Implementation of std::exception
	//@{

	/**	\brief Returns the cause of the exception (message)
	 **/
	virtual const char *what() const throw();

	//@}

private:
	void init(const char *clazz, const char *method,
	        const char *file, unsigned int line, const char *message) throw();
};


} // namespace libwfa

#endif // LIBWFA_LIBWFA_EXCEPTION_H

