#ifndef LIBWFA_INTERFACE_H
#define LIBWFA_INTERFACE_H

#include <stddef.h>

//! \name C-style interface to libwfa
//@{

/** \brief Perform analysis of state and difference density matrix (restricted)
    \param name Name of state (useable as filename)
    \param desc Description of state (one-line comment)
    \param ddm Difference density matrix in AO
    \param dm0 Ground state density matrix in AO
    \param nao Number of AO basis functions
 **/
void analyse_opdm(const char *name, const char *desc,
    double *ddm, double *dm0, size_t nao);

/** \brief Perform analysis of state and difference density matrix
        (unrestricted)
    \param name Name of state (useable as filename)
    \param desc Description of state (one-line comment)
    \param ddm_a Difference density matrix in AO (alpha spin)
    \param ddm_b Difference density matrix in AO (beta spin)
    \param dm0_a Ground state density matrix in AO (alpha spin)
    \param dm0_b Ground state density matrix in AO (beta spin)
    \param nao Number of AO basis functions
 **/
void analyse_opdm(const char *name, const char *desc,
    double *ddm_a, double *ddm_b, double *dm0_a, double *dm0_b, size_t nao);

/** \brief Perform analysis of state density matrix (restricted)
    \param name Name of state (useable as filename)
    \param desc Description of state (one-line comment)
    \param dm State density matrix in AO
    \param nao Number of AO basis functions
 **/
void analyse_opsdm(const char *name, const char *desc, double *dm, size_t nao);

/** \brief Perform analysis of state density matrix (unrestricted)
    \param name Name of state (useable as filename)
    \param desc Description of state (one-line comment)
    \param dm_a State density matrix in AO (alpha spin)
    \param dm_b State density matrix in AO (beta spin)
    \param nao Number of AO basis functions
 **/
void analyse_opsdm(const char *name, const char *desc,
    double *dm_a, double *dm_b, size_t nao);

/** \brief Perform analysis of transition density matrix (restricted)
    \param name Name of state (useable as filename)
    \param desc Description of state (one-line comment)
    \param tdm Transition density matrix in AO
    \param nao Number of AO basis functions
 **/
void analyse_optdm(const char *name, const char *desc, double *tdm, size_t nao);

/** \brief Perform analysis of transition density matrix (unrestricted)
    \param name Name of state (useable as filename)
    \param desc Description of state (one-line comment)
    \param tdm_a Transition density matrix in AO (alpha spin)
    \param tdm_b Transition density matrix in AO (beta spin)
    \param nao Number of AO basis functions
 **/
void analyse_optdm(const char *name, const char *desc,
    double *tdm_a, double *tdm_b, size_t nao);

/** \brief Perform SA-NTO analysis of transition density matrix (restricted)
    \param tdm Transition density matrix in AO
    \param nao Number of AO basis functions
 **/
void post_process_optdm(double *tdm, size_t nao);

/** \brief Perform SA-NTO analysis of transition density matrix (unrestricted)
    \param tdm_a Transition density matrix in AO (alpha spin)
    \param tdm_b Transition density matrix in AO (beta spin)
    \param nao Number of AO basis functions
 **/
void post_process_optdm(double *tdm_a, double *tdm_b, size_t nao);

/** \brief Reset the analysis to its original state after initialization
 **/
void reset_wf_analysis();

/** \brief Shutdown the analysis
 **/
void shutdown_wf_analysis();

//@}

#endif // LIBWFA_INTERFACE_H
