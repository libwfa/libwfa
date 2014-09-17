#ifndef LIBWFA_LIBWFA_H
#define LIBWFA_LIBWFA_H

/** \mainpage Wave-function analysis tool library

    libwfa provides a C++ library for the analysis of wavefunctions. This
    analysis is based on density matrices, which provide a unified formalism
    independent of the underlying wavefunction model.

    The methods used have been described in detail in Ref. [1] and examples
    are provided in Ref. [2].
    Briefly, the implemented methods consist of the following parts:

    1. Analysis of state density matrices for
       - population analysis,
       - density plotting,
       - natural orbitals, and
       - analysis of unpaired electrons.

    2. Analysis of transition density matrices for
       - plotting of the \a hole and \a electron densities,
       - plotting the transition density,
       - natural transition orbitals (with a possibility for state-averaging)
       - \a electron-hole correlation analysis using the charge transfer
           numbers (see also Ref. [3]).

    3. Analysis of difference density matrices for
       - plotting of the attachment/detachment densities,
       - attachment/detachment population analysis, and
       - natural difference orbitals.

    The required input from the quantum chemical program consists of the
    density matrices of interest, the AO-overlap matrix, the MO coefficients
    (used for orthogonalization), and information for population analysis.

    Current developments: multipole analysis.

    \b Literature

    [1] Plasser, F.; Wormit, M.; Dreuw, A. \a JCP, \a 2014, \a 141, 024106
        (DOI: 10.1063/1.4885819)
    [2] Plasser, F.; Baeppler, S.A.; Wormit, M.; Dreuw, A. \a JCP, \a 2014,
        \a 141, 024107 (DOI: 10.1063/1.4885820)
    [3] Plasser, F.; Lischka, H. \a JCTC, \b 2012, \a 8, 2777.
 **/

/** \defgroup libwfa Library of wave-function analysis tools
 **/

//! \name General
//@{
#include "version.h"
#include "wf_analysis.h"
//@}

//! \name Core interfaces
//@{
#include "core/mom_builder.h"
//@}

//! \name Various density matrix analyses
//@{
#include "analyses/ctnum_analysis.h"
#include "analyses/pop_mulliken.h"
#include "analyses/pop_loewdin.h"
//@}

//! \name Export and printing
//@{
#include "export/ctnum_export.h"
#include "export/density_printer_basic.h"
#include "export/density_printer_cube.h"
#include "export/density_printer_nil.h"
#include "export/export_cube_base.h"
#include "export/export_cube_i.h"
#include "export/export_molden_i.h"
#include "export/orbital_printer_basic.h"
#include "export/orbital_printer_cube.h"
#include "export/orbital_printer_molden.h"
#include "export/orbital_printer_nil.h"
//@}

#endif // LIBWFA_LIBWFA_H
