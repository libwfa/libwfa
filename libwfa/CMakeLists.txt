set(HEADERS
    ab_matrix.h
    ab_vector.h
    analyse_sdm.h
    analyse_tdm.h
    ctnum_analysis_i.h
    ctnum_analysis.h
    ctnum_export.h
    ctnum_data_i.h
    ctnumbers.h
    density_type.h
    ev_data_i.h
    export_cube_i.h
    export_molden_i.h
    export_densities_cube.h
    export_densities_i.h
    export_orbitals_i.h
    export_orbitals_cube.h
    export_orbitals_molden.h
    libwfa_exception.h
    ndo_analysis.h
    ndo_data_print.h
    nto_analysis.h
    nto_data_print.h
    no_analysis.h
    no_data_print.h
    orbital_type.h
    pop_analysis_ad.h
    pop_analysis_dm.h
    pop_analysis_i.h
    pop_data.h
    pop_mulliken.h
    pop_print_default.h
    pop_print_i.h
    transformations_dm.h
    selector.h
    state_info.h
    version.h
)

set(SRC
    analyse_sdm.C
    analyse_tdm.C
    ctnum_analysis.C
    ctnum_export.C
    ctnumbers.C
    density_type.C
    export_densities_cube.C
    export_orbitals_cube.C
    export_orbitals_molden.C
    libwfa_exception.C
    ndo_analysis.C
    ndo_data_print.C
    nto_analysis.C
    nto_data_print.C
    no_analysis.C
    no_data_print.C
    orbital_type.C
    pop_analysis_ad.C
    pop_analysis_dm.C
    pop_mulliken.C
    pop_print_default.C
    transformations_dm.C
    selector.C
    state_info.C
    version.C
)

set(HEADERS_QCHEM
    qchem/export_cube_qchem.h
    qchem/export_molden_qchem.h
)

set(SRC_QCHEM
    qchem/export_cube_qchem.C
    qchem/export_molden_qchem.C
)


if(QCHEM_SOURCE_DIR)
    add_library(wfa_qchem ${SRC} ${HEADERS} ${SRC_QCHEM} ${HEADERS_QCHEM})
    set(TARGETS wfa_qchem)
else(QCHEM_SOURCE_DIR)
    add_library(wfa ${SRC} ${HEADERS})
    set(TARGETS wfa)
endif(QCHEM_SOURCE_DIR)


