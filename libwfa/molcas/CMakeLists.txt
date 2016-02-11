set(HDR
    molcas_wf_analysis_data.h
)

set(SRC
    molcas_wf_analysis_data.C
    molcas_state_analysis.C
)

set(TARGETS wfa_molcas.x)
add_executable(wfa_molcas.x ${SRC} ${HDR})
add_dependencies(wfa_molcas.x generated_headers)
target_link_libraries(wfa_molcas.x wfa)
target_link_libraries(wfa_molcas.x hdf5 hdf5_cpp)

set_property(TARGET ${TARGETS} APPEND PROPERTY 
    COMPILE_DEFINITIONS_DEBUG LIBWFA_DEBUG)