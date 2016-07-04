
libxfts_sm_a_SOURCES = $(addprefix $(topsrcdir)/sm/, \
	xfts_sparse_matrix_mod.F90 \
	xfts_sparse_matrix_type.F90 \
	fts_qsortc.c \
	fts_sparse_matrix_enum.F90  )



libcfts_sm_a_SOURCES = $(subst xfts_,cfts_,$(libxfts_sm_a_SOURCES))
libdfts_sm_a_SOURCES = $(subst xfts_,dfts_,$(libxfts_sm_a_SOURCES))
libsfts_sm_a_SOURCES = $(subst xfts_,sfts_,$(libxfts_sm_a_SOURCES))
libzfts_sm_a_SOURCES = $(subst xfts_,zfts_,$(libxfts_sm_a_SOURCES))

libcfts_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcfts_sm_a_SOURCES)))
libdfts_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdfts_sm_a_SOURCES)))
libsfts_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsfts_sm_a_SOURCES)))
libzfts_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzfts_sm_a_SOURCES)))

libcfts_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libcfts_sm_a_SOURCES)))
libdfts_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libdfts_sm_a_SOURCES)))
libsfts_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libsfts_sm_a_SOURCES)))
libzfts_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libzfts_sm_a_SOURCES)))

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcfts_sm_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdfts_sm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsfts_sm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzfts_sm_a_OBJECTS:.o=.d)
endif

