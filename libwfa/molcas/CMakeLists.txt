set(HDR
    molcas_wf_analysis.h
    molcas_wf_analysis_data.h
    molcas_mom_builder.h
    molcas_export_h5orbs.h
)

set(SRC
    molcas_wf_analysis.C
    molcas_wf_analysis_data.C
    molcas_state_analysis.C
    molcas_mom_builder.C
    molcas_export_h5orbs.C
)

set(TARGETS wfa_molcas.x)
add_executable(wfa_molcas.x ${SRC} ${HDR})
# add_dependencies(wfa_molcas.x generated_headers)

set_property(TARGET ${TARGETS} APPEND PROPERTY 
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)
    
target_link_libraries(wfa_molcas.x wfa ${LIBWFA_DEP_LIBS})
target_link_libraries(wfa_molcas.x hdf5 hdf5_cpp)
