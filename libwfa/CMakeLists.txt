set(ANALYSES_HEADERS
    analyses/ctnum_analysis_i.h
    analyses/ctnum_analysis.h
    analyses/ctnumbers.h
    analyses/exciton_analysis_ad.h
    analyses/exciton_analysis_base.h
    analyses/exciton_analysis.h
    analyses/exciton_moments.h
    analyses/ndo_analysis.h
    analyses/nto_analysis.h
    analyses/no_analysis.h
    analyses/pop_analysis_ad.h
    analyses/pop_analysis_dm.h
    analyses/pop_analysis_i.h
    analyses/pop_data.h
    analyses/pop_loewdin.h
    analyses/pop_mulliken.h
    analyses/sa_nto_analysis.h
)

set(ANALYSES_SRC
    analyses/ctnum_analysis.C
    analyses/ctnumbers.C
    analyses/exciton_analysis_ad.C
    analyses/exciton_analysis_base.C
    analyses/exciton_analysis.C
    analyses/ndo_analysis.C
    analyses/nto_analysis.C
    analyses/no_analysis.C
    analyses/pop_analysis_ad.C
    analyses/pop_analysis_dm.C
    analyses/pop_data.C
    analyses/pop_loewdin.C
    analyses/pop_mulliken.C
    analyses/sa_nto_analysis.C
)


set(CORE_HEADERS
    core/ab_matrix.h
    core/ab_object.h
    core/ab_vector.h
    core/constants.h
    core/density_type.h
    core/mom_builder_i.h
    core/mom_builder.h
    core/orbital_data.h
    core/orbital_selector.h
    core/orbital_type.h
    core/selector.h
)

set(CORE_SRC
    core/constants.C
    core/density_type.C
    core/mom_builder.C
    core/orbital_data.C
    core/orbital_selector.C
    core/orbital_type.C
    core/selector.C
)


set(EXPORT_HEADERS
    export/ctnum_export.h
    export/ctnum_printer_i.h
    export/cube_writer.h
    export/density_printer_cube.h
    export/density_printer_i.h
    export/density_printer_none.h
    export/density_printer_print.h
    export/export_cube_base.h
    export/export_cube_i.h
    export/export_molden_i.h
    export/grid3d.h
    export/orbital_printer_basic.h
    export/orbital_printer_cube.h
    export/orbital_printer_i.h
    export/orbital_printer_molden.h
    export/orbital_printer_nil.h
    export/pop_printer_default.h
    export/pop_printer_i.h
)

set(EXPORT_SRC
    export/ctnum_export.C
    export/cube_writer.C
    export/density_printer_cube.C
    export/density_printer_print.C
    export/export_cube_base.C
    export/grid3d.C
    export/orbital_printer_basic.C
    export/orbital_printer_cube.C
    export/orbital_printer_molden.C
    export/pop_printer_default.C
)


set(HEADERS
    ${CORE_HEADERS}
    ${EXPORT_HEADERS}
    ${ANALYSES_HEADERS}
#    analyse_opdm.h
#    analyse_optdm.h
    libwfa_exception.h
    version.h
)

set(SRC
    ${CORE_SRC}
    ${EXPORT_SRC}
    ${ANALYSES_SRC}
#    analyse_opdm.C
#    analyse_optdm.C
    libwfa_exception.C
    version.C
)

set(TARGETS wfa)
add_library(wfa ${SRC} ${HEADERS})

set_property(TARGET ${TARGETS} APPEND PROPERTY 
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)

