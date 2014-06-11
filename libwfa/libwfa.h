#ifndef LIBWFA_LIBWFA_H
#define LIBWFA_LIBWFA_H

/** \mainpage Wave-function analysis tool library

    The library consists of a number of tools for wave function or rather
    density matrix analysis.

    ... to be continued ...
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
