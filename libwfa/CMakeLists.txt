set(ANALYSES_HEADERS
    analyses/ctnum_analysis_i.h
    analyses/ctnum_analysis.h
    analyses/ctnumbers.h
    analyses/ndo_analysis.h
    analyses/nto_analysis.h
    analyses/no_analysis.h
    analyses/pop_analysis_ad.h
    analyses/pop_analysis_dm.h
    analyses/pop_analysis_i.h
    analyses/pop_data.h
    analyses/pop_mulliken.h
)

set(ANALYSES_SRC
    analyses/ctnum_analysis.C
    analyses/ctnumbers.C
    analyses/ndo_analysis.C
    analyses/nto_analysis.C
    analyses/no_analysis.C
    analyses/pop_analysis_ad.C
    analyses/pop_analysis_dm.C
    analyses/pop_mulliken.C
)


set(CORE_HEADERS
    core/ab_matrix.h
    core/ab_object.h
    core/ab_selector.h
    core/ab_vector.h
    core/density_type.h
    core/orbital_type.h
    core/selector.h
    core/transformations_dm.h
)

set(CORE_SRC
    core/density_type.C
    core/orbital_type.C
    core/selector.C
    core/transformations_dm.C
)


set(EXPORT_HEADERS
    export/ctnum_export.h
    export/ctnum_printer_i.h
    export/cube_writer.h
    export/ev_printer_i.h
    export/export_densities_i.h
    export/export_densities_none.h
    export/export_densities_print.h
    export/export_orbitals_i.h
    export/export_orbitals_none.h
    export/export_orbitals_print.h
    export/grid3d.h
    export/ndo_data_print.h
    export/no_data_print.h
    export/nto_data_print.h
    export/pop_printer_default.h
    export/pop_printer_i.h
)

set(EXPORT_SRC
    export/ctnum_export.C
    export/cube_writer.C
    export/export_densities_print.C
    export/export_orbitals_print.C
    export/grid3d.C
    export/ndo_data_print.C
    export/nto_data_print.C
    export/no_data_print.C
    export/pop_printer_default.C
)


set(HEADERS
    ${CORE_HEADERS}
    ${EXPORT_HEADERS}
    ${ANALYSES_HEADERS}
    analyse_sdm.h
    analyse_tdm.h
    libwfa_exception.h
    version.h
)

set(SRC
    ${CORE_SRC}
    ${EXPORT_SRC}
    ${ANALYSES_SRC}
    analyse_sdm.C
    analyse_tdm.C
    libwfa_exception.C
    version.C
)

set(HEADERS_QCHEM
#    qchem/export_cube_qchem.h
#    qchem/export_molden_qchem.h
)

set(SRC_QCHEM
#    qchem/export_cube_qchem.C
#    qchem/export_molden_qchem.C
)

set(TARGETS wfa)
add_library(wfa ${SRC} ${HEADERS})

if(QCHEM_SOURCE_DIR)
    set(TARGETS ${TARGETS} wfa_qchem)
    add_library(wfa_qchem ${SRC_QCHEM} ${HEADERS_QCHEM})
endif(QCHEM_SOURCE_DIR)

set_property(TARGET ${TARGETS} APPEND PROPERTY 
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)

