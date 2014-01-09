set(SRC_HEADER
	libwfa_exception.h
	pop_analysis_i.h
	pop_data.h
	pop_mulliken.h
	pop_print_default.h
	pop_print_i.h
    version.h
)

set(SRC
	libwfa_exception.C
	pop_mulliken.C
	pop_print_default.C
    version.C
)

add_library(wfa ${SRC} ${SRC_HEADER})


