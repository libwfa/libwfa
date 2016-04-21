set(HDR
    molcas_wf_analysis.h
    molcas_wf_analysis_data.h
    molcas_mom_builder.h
    molcas_export_h5orbs.h
)

set(SRC
    main.C
    molcas_wf_analysis.C
    molcas_wf_analysis_data.C
    molcas_mom_builder.C
    molcas_export_h5orbs.C
)

set(TARGETS wfa.exe)
add_executable(wfa.exe ${SRC} ${HDR})
# add_dependencies(wfa_molcas.x generated_headers)

set_property(TARGET ${TARGETS} APPEND PROPERTY 
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)
    
target_link_libraries(wfa.exe wfa ${LIBWFA_DEP_LIBS})
target_link_libraries(wfa.exe hdf5 hdf5_cpp)
