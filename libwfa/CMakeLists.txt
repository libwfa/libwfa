set(SRC_HEADER
	ab_matrix.h
	ab_vector.h
	ctnum_analysis_i.h
	ctnum_analysis.h
	ctnum_data.h
	ctnum_export.h
	ctnum_print_i.h
	libwfa_exception.h
	pop_analysis_i.h
	pop_data.h
	pop_mulliken.h
	pop_print_default.h
	pop_print_i.h
    transformations_dm.h
    selector.h
    version.h
)

set(SRC
	ctnum_analysis.C
	ctnum_export.C
	libwfa_exception.C
	pop_mulliken.C
	pop_print_default.C
    transformations_dm.C
    selector.C
    version.C
)

add_library(wfa ${SRC} ${SRC_HEADER})


