#ifndef LIBWFA_LIBWFA_H
#define LIBWFA_LIBWFA_H

/** \mainpage Wave-function analysis tool library

    The library consists of a number of tools for wave function analysis.

    ... to be continued ...
 **/

/** \defgroup libwfa Library of wave-function analysis tools
 **/

#include "version.h"
#include "analyse_opdm.h"
#include "analyse_optdm.h"

#include "ctnumbers.h"
#include "ndo_analysis.h"
#include "no_analysis.h"
#include "nto_analysis.h"
#include "pop_analysis_ad.h"
#include "pop_analysis_dm.h"

#include "ctnum_analysis.h"
#include "pop_mulliken.h"

#include "ctnum_export.h"
#include "pop_printer_default.h"

#include "ev_printer_no.h"
#include "ev_printer_none.h"
#include "ev_printer_ndo.h"
#include "ev_printer_nto.h"

#include "export_data_cm.h"
#include "export_data_cube.h"
#include "export_data_none.h"
#include "export_data_print.h"
#include "export_orbitals_molden.h"

#endif // LIBWFA_LIBWFA_H
