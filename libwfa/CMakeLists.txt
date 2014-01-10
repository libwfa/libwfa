set(SRC_HEADER
	ab_matrix.h
	ab_vector.h
	libwfa_exception.h
	pop_analysis_i.h
	pop_data.h
	pop_mulliken.h
	pop_print_default.h
	pop_print_i.h
    transformations_dm.h
    version.h
)

set(SRC
	libwfa_exception.C
	pop_mulliken.C
	pop_print_default.C
    transformations_dm.C
    version.C
)

add_library(wfa ${SRC} ${SRC_HEADER})


