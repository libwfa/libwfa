#ifndef LIBWFA_LIBWFA_H
#define LIBWFA_LIBWFA_H

/** \mainpage Wave-function analysis tool library

    libwfa provides a C++ library for the analysis of wavefunctions. This
    analysis is based on density matrices, which provide a unified formalism
    independent of the wavefunction model.

    The methods used have been described in detail in Ref. [1] and examples are
    provided in Ref. [2].
    Briefly, the implemented methods consist of the following parts:

    1. Analysis of state density matrices for
       - population analysis,
       - density plotting,
       - natural orbitals, and
       - analysis of unpaired electrons.

    2. Analysis of transition density matrices for
       - plotting of the \a hole and \a electron densities,
       - plotting the transition density,
       - natural transition orbitals (with a possibility for state-averaging), and
       - \a electron-hole correlation analysis using the charge transfer numbers (see also Ref. [3]).

    3. Analysis of difference density matrices for
       - plotting of the attachment/detachment densities,
       - attachment/detachment population analysis, and
       - natural difference orbitals.

    The required input from the quantum chemical program consists of the
    density matrices of interest, the AO-overlap matrix, the MO coefficients
    (used for orthogonalization), and information for population analysis.

    Current developments: multipole analysis.

    \b Literature

    [1] Plasser, F.; Wormit, M.; Dreuw, A. <em>submitted for publication</em>

    [2] Plasser, F.; Baeppler, S.A.; Wormit, M.; Dreuw, A. <em>submitted for publication</em>

    [3] Plasser, F.; Lischka, H. \a JCTC, \b 2012, \a 8, 2777.
 **/

/** \defgroup libwfa Library of wave-function analysis tools
 **/

//! \name Various density matrix analyses
//@{
#include "analyses/ctnum_analysis.h"
#include "analyses/ctnumbers.h"
#include "analyses/ndo_analysis.h"
#include "analyses/no_analysis.h"
#include "analyses/nto_analysis.h"
#include "analyses/pop_analysis_ad.h"
#include "analyses/pop_analysis_dm.h"
#include "analyses/pop_mulliken.h"
//@}

//! \name Export and printing
//@{
#include "export/ctnum_export.h"
#include "export/ev_printer_no.h"
#include "export/ev_printer_none.h"
#include "export/ev_printer_ndo.h"
#include "export/ev_printer_nto.h"
#include "export/export_data_cm.h"
#include "export/export_data_cube.h"
#include "export/export_data_none.h"
#include "export/export_data_print.h"
#include "export/export_orbitals_molden.h"
#include "export/pop_printer_default.h"
//@}

//! \name General
//@{
#include "version.h"
#include "analyse_opdm.h"
#include "analyse_optdm.h"
//@}


#endif // LIBWFA_LIBWFA_H
