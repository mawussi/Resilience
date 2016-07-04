

libxfts_dm_a_SOURCES = $(topsrcdir)/dm/xfts_dense_matrix_mod.F90 $(topsrcdir)/dm/xfts_dense_matrix_type.F90

libcfts_dm_a_SOURCES = $(subst xfts_,cfts_,$(libxfts_dm_a_SOURCES))
libdfts_dm_a_SOURCES = $(subst xfts_,dfts_,$(libxfts_dm_a_SOURCES))
libsfts_dm_a_SOURCES = $(subst xfts_,sfts_,$(libxfts_dm_a_SOURCES))
libzfts_dm_a_SOURCES = $(subst xfts_,zfts_,$(libxfts_dm_a_SOURCES))

libfts_dm_a_SOURCES = $(topsrcdir)/dm/fts_dense_matrix_mod.F90                 \
	$(sort $(libcfts_dm_a_SOURCES) $(libdfts_dm_a_SOURCES) \
	$(libsfts_dm_a_SOURCES) $(libzfts_dm_a_SOURCES))

libcfts_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcfts_dm_a_SOURCES)))
libdfts_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdfts_dm_a_SOURCES)))
libsfts_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsfts_dm_a_SOURCES)))
libzfts_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzfts_dm_a_SOURCES)))

libfts_dm_a_OBJECTS  += $(patsubst %.F90,%.o,$(filter %.F90,$(libfts_dm_a_SOURCES))) 

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcfts_dm_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdfts_dm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsfts_dm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzfts_dm_a_OBJECTS:.o=.d)
endif


