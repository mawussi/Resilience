
libxfts_part_a_SOURCES = $(addprefix $(topsrcdir)/part/, \
	fts_part_type.F90 \
	fts_get_neighbors_mod.F90 \
	 xfts_analyze_mod.F90 )


	# mph_domain_mod.F90 \
	# mph_rhs_partition_mod.F90 \
	# xmph_part_mod.F90 \
	# xmph_part_ordermatrix_mod.F90 \
	# xmph_part_builddomains_mod.F90 \
	# xmph_part_distmatrix_mod.F90 \
	# xmph_part_distrhs_mod.F90 \
	# xmph_part_collectsol_mod.F90 )

# Special case for METIS.
# Note: this is uneccessary in for mph_part_scotch_mod.F90 which already handles internally -DHAVE_LIBSCOTCH
# ifneq ($(filter $(METIS_FCFLAGS),-DHAVE_LIBMETIS),)
# libxmph_part_a_SOURCES += $(addprefix $(topsrcdir)/part/, \
# 	mph_part_metis.c mph_part_metis_fortran.c \
# 	mph_part_metis_proto.h mph_part_metis_struct.h )
# endif


libcfts_part_a_SOURCES = $(subst xfts_,cfts_,$(libxfts_part_a_SOURCES))
libdfts_part_a_SOURCES = $(subst xfts_,dfts_,$(libxfts_part_a_SOURCES))
libsfts_part_a_SOURCES = $(subst xfts_,sfts_,$(libxfts_part_a_SOURCES))
libzfts_part_a_SOURCES = $(subst xfts_,zfts_,$(libxfts_part_a_SOURCES))

libcfts_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcfts_part_a_SOURCES)))
libdfts_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdfts_part_a_SOURCES)))
libsfts_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsfts_part_a_SOURCES)))
libzfts_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzfts_part_a_SOURCES)))

libcfts_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libcfts_part_a_SOURCES)))
libdfts_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libdfts_part_a_SOURCES)))
libsfts_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libsfts_part_a_SOURCES)))
libzfts_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libzfts_part_a_SOURCES)))

# Warning : order in EXTRA_*CFLAGS defined below are important.
#
# Since we call specific subfunctions of libmetis.a.
# we must include metis.h from METIS before the one defined by SCOTCH.
#
#
$(libcfts_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(libdfts_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(libsfts_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(libzfts_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 

$(libcfts_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)
$(libdfts_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)
$(libsfts_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)
$(libzfts_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcfts_part_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdfts_part_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsfts_part_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzfts_part_a_OBJECTS:.o=.d)
endif

