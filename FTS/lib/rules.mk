# @file Makefile
#
# MaPHyS is a software package provided by INRIA.
#
# Build libmaphys.a or lib[cdsz]maphys.a, and copy its dependencies
#
# @author Yohan Lee-tin-yien
#
###

libcftsolver_a_LIBADD = $(libfts_common_a_OBJECTS) \
	$(libcfts_dm_a_OBJECTS) $(libcfts_sm_a_OBJECTS) \
	$(libcfts_part_a_OBJECTS) $(libcfts_src_a_OBJECTS) 


libdftsolver_a_LIBADD = $(libfts_common_a_OBJECTS) \
	$(libdfts_dm_a_OBJECTS) $(libdfts_sm_a_OBJECTS) \
	$(libdfts_part_a_OBJECTS) $(libdfts_src_a_OBJECTS) 


libsftsolver_a_LIBADD = $(libfts_common_a_OBJECTS) \
	$(libsfts_dm_a_OBJECTS) $(libsfts_sm_a_OBJECTS) \
	$(libsfts_part_a_OBJECTS) $(libsfts_src_a_OBJECTS) 


libzftsolver_a_LIBADD = $(libfts_common_a_OBJECTS) \
	$(libzfts_dm_a_OBJECTS) $(libzfts_sm_a_OBJECTS) \
	$(libzfts_part_a_OBJECTS) $(libzfts_src_a_OBJECTS) 

##
libftsolver_a_LIBADD = 

ifneq (,$(filter c,$(ARITHS)))
libftsolver_a_LIBADD += $(libcftsolver_a_LIBADD)
endif

ifneq (,$(filter d,$(ARITHS)))
libftsolver_a_LIBADD += $(libdftsolver_a_LIBADD)
endif

ifneq (,$(filter z,$(ARITHS)))
libftsolver_a_LIBADD += $(libsftsolver_a_LIBADD)
endif

ifneq (,$(filter z,$(ARITHS)))
libftsolver_a_LIBADD += $(libzftsolver_a_LIBADD)
endif
