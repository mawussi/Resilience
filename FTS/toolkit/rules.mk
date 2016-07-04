

libxfts_toolkit_a_SOURCES = $(addprefix $(topsrcdir)/toolkit/, \
	fts_toolkit_bindproc.c \
	fts_toolkit_parser_mod.F90 \
	xfts_toolkit_read_mod.F90 \
	xfts_toolkit_gen_mod.F90 \
	xfts_toolkit_mod.F90 )


libcfts_toolkit_a_LIBADD =$(libfts_common_a_OBJECTS)  $(libcfts_dm_a_OBJECTS) $(libcfts_sm_a_OBJECTS) 
libdfts_toolkit_a_LIBADD =$(libfts_common_a_OBJECTS) $(libdfts_dm_a_OBJECTS) $(libdfts_sm_a_OBJECTS) 
libsfts_toolkit_a_LIBADD =$(libfts_common_a_OBJECTS) $(libsfts_dm_a_OBJECTS) $(libsfts_sm_a_OBJECTS) 
libzfts_toolkit_a_LIBADD =$(libfts_common_a_OBJECTS) $(libzfts_dm_a_OBJECTS) $(libzfts_sm_a_OBJECTS) 

### 
libcfts_toolkit_a_SOURCES = $(subst xfts_,cfts_,$(libxfts_toolkit_a_SOURCES))
libdfts_toolkit_a_SOURCES = $(subst xfts_,dfts_,$(libxfts_toolkit_a_SOURCES))
libsfts_toolkit_a_SOURCES = $(subst xfts_,sfts_,$(libxfts_toolkit_a_SOURCES))
libzfts_toolkit_a_SOURCES = $(subst xfts_,zfts_,$(libxfts_toolkit_a_SOURCES))

libcfts_toolkit_a_OBJECTS = \
	$(patsubst %.F90,%.o,$(filter %.F90,$(libcfts_toolkit_a_SOURCES))) \
	$(patsubst %.c,%.o,$(filter %.c,$(libcfts_toolkit_a_SOURCES)))

libdfts_toolkit_a_OBJECTS = \
	$(patsubst %.F90,%.o,$(filter %.F90,$(libdfts_toolkit_a_SOURCES))) \
	$(patsubst %.c,%.o,$(filter %.c,$(libdfts_toolkit_a_SOURCES)))

libsfts_toolkit_a_OBJECTS = \
	$(patsubst %.F90,%.o,$(filter %.F90,$(libsfts_toolkit_a_SOURCES))) \
	$(patsubst %.c,%.o,$(filter %.c  ,$(libsfts_toolkit_a_SOURCES)))

libzfts_toolkit_a_OBJECTS = \
	$(patsubst %.F90,%.o,$(filter %.F90,$(libzfts_toolkit_a_SOURCES))) \
	$(patsubst %.c,%.o,$(filter %.c,$(libzfts_toolkit_a_SOURCES)))



# $(libcfts_toolkit_a_OBJECTS) : EXTRA_CFLAGS=$(HWLOC_CFLAGS)
# $(libdfts_toolkit_a_OBJECTS) : EXTRA_CFLAGS=$(HWLOC_CFLAGS)
# $(libsfts_toolkit_a_OBJECTS) : EXTRA_CFLAGS=$(HWLOC_CFLAGS)
# $(libzfts_toolkit_a_OBJECTS) : EXTRA_CFLAGS=$(HWLOC_CFLAGS)

# $(libcfts_toolkit_a_OBJECTS) : | $(libcfts_toolkit_a_DEPENDENCIES)
# $(libdfts_toolkit_a_OBJECTS) : | $(libdfts_toolkit_a_DEPENDENCIES)
# $(libsfts_toolkit_a_OBJECTS) : | $(libsfts_toolkit_a_DEPENDENCIES)
# $(libzfts_toolkit_a_OBJECTS) : | $(libzfts_toolkit_a_DEPENDENCIES)

#
ifneq (,$(filter c,$(ARITHS)))
include $(libcfts_toolkit_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter d,$(ARITHS)))
include $(libdfts_toolkit_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsfts_toolkit_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzfts_toolkit_a_OBJECTS:.o=.d)
endif
