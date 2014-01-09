set(SRC_HEADER
	pop_analysis_i.h
	pop_mulliken.h
	pop_print_default.h
	pop_print_i.h
)

set(SRC
	pop_mulliken.C
	pop_print_default.C
)

add_library(wfa ${SRC} ${SRC_HEADER})


