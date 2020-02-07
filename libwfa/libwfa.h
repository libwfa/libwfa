#ifndef LIBWFA_LIBWFA_H
#define LIBWFA_LIBWFA_H

/** \mainpage Wave-function analysis tool library

    \c libwfa provides a C++ library for the analysis of wavefunctions. This
    analysis is based on density matrices, which provide a unified formalism
    independent of the underlying wavefunction model.

    The methods which have been implemented are described in detail
    in Ref. [1] and examples are provided in Ref. [2].
    Further developments are covered in Refs [4] and [5].
    Briefly, the implemented methods consist of the following parts:

    1. Analysis of state density matrices for
       - population analysis,
       - density plotting,
       - natural orbitals, and
       - analysis of unpaired electrons.

    2. Analysis of transition density matrices for
       - plotting of the \a hole and \a electron densities,
       - plotting the transition density,
       - natural transition orbitals (NTOs) with a possibility for state-averaging,
       - entanglement analysis based on NTOs [6],
       - \a electron-hole correlation analysis using the charge transfer
           numbers (see also Ref. [3]).

    3. Analysis of difference density matrices for
       - plotting of the attachment/detachment densities,
       - attachment/detachment population analysis, and
       - natural difference orbitals.

    4. Real-space analysis of properties within an exciton picture,
       see Refs [4], [5].

    5. Electrostatic potentials of various effective densities.

    The required input from the quantum chemical program consists of the
    density matrices of interest, the AO-overlap matrix, the MO coefficients
    (used for orthogonalization), and information for population analysis.
    For the exciton analysis routines also the dipole and quadrupole integrals
    are needed.

    \par Literature
    -# Plasser, F.; Wormit, M.; Dreuw, A. \a JCP, \b 2014, \a 141, 024106
        (DOI: 10.1063/1.4885819).
    -# Plasser, F.; Baeppler, S.A.; Wormit, M.; Dreuw, A. \a JCP, \b 2014,
        \a 141, 024107 (DOI: 10.1063/1.4885820).
    -# Plasser, F.; Lischka, H. \a JCTC, \b 2012, \a 8, 2777.
    -# S. Baeppler, F. Plasser, M. Wormit, A. Dreuw \a PRA, \b 2014, \a 90, 052521
        (DOI: 10.1103/PhysRevA.90.052521).
    -# F. Plasser, B. Thomitzni, S. Baeppler, J. Wenzel, D. Rehn, M. Wormit, A. Dreuw
        \a JCC, \b 2015, \a 36, 1609-1620 (DOI: 10.1002/jcc.23975).
    -# F. Plasser \a JCP, \b 2016, \a 144, 194107 (DOI: 10.1063/1.4949535).

    \author Felix Plasser
    \author Michael Wormit
    \author Stefanie Mewes (Baeppler)
    \author Benjamin Thomitzni
    \author Feng Chen
    \author Anna I. Krylov
    \author Pavel Pokhilko

    \date 2014-2019

    \copyright (c) 2014, F. Plasser and M. Wormit
    \copyright All rights reserved.
    \copyright Modifications copyright (C) 2019, Loughborough University.
    \copyright Redistribution and use in source and binary forms, with or
        without modification, are permitted provided that the following
        conditions are met:
    \copyright
     1. Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
     3. Neither the name of the copyright holder nor the names of its
        contributors may be used to endorse or promote products derived from
        this software without specific prior written permission.
    \copyright
        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
        "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
        LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
        A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
        HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
        LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
        DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
        THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
        OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **/

/** \defgroup libwfa Library of wave-function analysis tools
 **/

//! \name General
//@{
#include "version.h"
#include "libwfa_exception.h"
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
