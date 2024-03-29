set(ANALYSES_HEADERS
    analyses/ctnum_analysis_i.h
    analyses/ctnum_analysis.h
    analyses/ctnumbers.h
    analyses/dyson_analysis.h
    analyses/dens_mom.h
    analyses/exciton_analysis_ad.h
    analyses/exciton_analysis_base.h
    analyses/exciton_analysis.h
    analyses/exciton_moments.h
    analyses/ndo_analysis.h
    analyses/nto_analysis.h
    analyses/no_analysis.h
    analyses/pop_analysis_ad.h
    analyses/pop_analysis_dm.h
    analyses/pop_analysis_tdm.h
    analyses/pop_analysis_i.h
    analyses/pop_data.h
    analyses/pop_loewdin.h
    analyses/pop_mulliken.h
    analyses/sa_nto_analysis.h
    analyses/om_descriptor.h
        )

set(ANALYSES_SRC
    analyses/ctnum_analysis.C
    analyses/ctnumbers.C
    analyses/dyson_analysis.C
    analyses/dens_mom.C
    analyses/exciton_analysis_ad.C
    analyses/exciton_analysis_base.C
    analyses/exciton_analysis.C
    analyses/ndo_analysis.C
    analyses/nto_analysis.C
    analyses/no_analysis.C
    analyses/pop_analysis_ad.C
    analyses/pop_analysis_dm.C
    analyses/pop_analysis_tdm.C
    analyses/pop_data.C
    analyses/pop_loewdin.C
    analyses/pop_mulliken.C
    analyses/sa_nto_analysis.C
    analyses/om_descriptor.C
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
    export/density_printer_basic.h
    export/density_printer_cube.h
    export/density_printer_i.h
    export/density_printer_nil.h
    export/export_cube_base.h
    export/export_cube_i.h
    export/export_molden_i.h
    export/grid3d.h
    export/orbital_printer_basic.h
    export/orbital_printer_cube.h
    export/orbital_printer_i.h
    export/orbital_printer_molden.h
    export/orbital_printer_nil.h
)

set(EXPORT_SRC
    export/ctnum_export.C
    export/cube_writer.C
    export/density_printer_basic.C
    export/density_printer_cube.C
    export/export_cube_base.C
    export/grid3d.C
    export/orbital_printer_basic.C
    export/orbital_printer_cube.C
    export/orbital_printer_molden.C
)

set(HEADERS
    ${CORE_HEADERS}
    ${EXPORT_HEADERS}
    ${ANALYSES_HEADERS}
    libwfa_exception.h
    libwfa_interface.h
    version.h
    wf_analysis_data_i.h
    wf_analysis.h
)

set(SRC
    ${CORE_SRC}
    ${EXPORT_SRC}
    ${ANALYSES_SRC}
    libwfa_exception.C
    libwfa_interface.C
    version.C
    wf_analysis.C
)

set(MOLCAS_HDR
    molcas/molcas_wf_analysis.h
    molcas/molcas_wf_analysis_data.h
    molcas/molcas_mom_builder.h
    molcas/molcas_export_h5orbs.h
)

set(MOLCAS_SRC
    molcas/driver.C
    molcas/molcas_wfa_interface.C
    molcas/molcas_wf_analysis.C
    molcas/molcas_wf_analysis_data.C
    molcas/molcas_mom_builder.C
    molcas/molcas_export_h5orbs.C
)

if (CMAKE_VERSION VERSION_LESS "3.1")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

if(MOLCAS_LIB)
  Message("Configuring libwfa with Molcas objects included ...")
  project(libwfa CXX)
  add_definitions ( -DARMA_BLAS_LONG )
  add_library(wfa_molcas ${SRC} ${HEADERS} ${MOLCAS_HDR} ${MOLCAS_SRC})
  install (TARGETS wfa_molcas
           RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
           LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
else()
    set(TARGETS wfa)
    add_library(wfa STATIC ${SRC} ${HEADERS})

    if(QCHEM_SOURCE_DIR)
      target_link_libraries(wfa armadillo)
      add_subdirectory(qchem)
    endif(QCHEM_SOURCE_DIR)

    if(MOLCAS_EXE)
      Message("Configuring standalone Molcas executable ...")
      add_subdirectory(molcas)
    endif(MOLCAS_EXE)

endif(MOLCAS_LIB)


set_property(TARGET ${TARGETS} APPEND PROPERTY
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)
