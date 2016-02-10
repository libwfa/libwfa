set(HDR
    molcas_wf_analysis_data.h
)

set(SRC
    molcas_wf_analysis_data.C
    molcas_state_analysis.C
)

set(TARGETS wfa_molcas)

add_executable(wfa_molcas.x ${SRC} ${HDR})
TARGET_LINK_LIBRARIES(wfa_molcas.x hdf5 hdf5_cpp)
