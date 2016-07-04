
libxfts_src_a_SOURCES = $(addprefix $(topsrcdir)/src/, \
	xfts_ftsolver_type.F90 xfts_ftsolver_mod.F90  \
	xfts_mumps_mod.F90 xfts_precond_mod.F90  \
	xfts_sls_mod.F90 fts_ftsolver_enum.F90  \
	xfts_compute_aux_matrices_mod.F90 xfts_solve_mod.F90 \
	xfts_fault_mod.F90 xfts_state_print_mod.F90 xfts_state_update_mod.F90 )

xfts_test_LDADD := $(addprefix $(buildlibdir)/,\
	libftsolver.a   libxfts_toolkit.a libslatec.a libpackgmres.a libpackfgmres.a)

xfts_bench_qrm_LDADD := $(addprefix $(buildlibdir)/,\
	libftsolver.a   libxfts_toolkit.a libslatec.a )



xfts_test_SOURCES :=  $(topsrcdir)/src/xfts_test.F90
xfts_bench_qrm_SOURCES :=  $(topsrcdir)/src/xfts_bench_qrm.F90


cfts_test_SOURCES := $(subst xfts_,cfts_,$(xfts_test_SOURCES))
dfts_test_SOURCES := $(subst xfts_,dfts_,$(xfts_test_SOURCES))
sfts_test_SOURCES := $(subst xfts_,sfts_,$(xfts_test_SOURCES))
zfts_test_SOURCES := $(subst xfts_,zfts_,$(xfts_test_SOURCES))

cfts_bench_qrm_SOURCES := $(subst xfts_,cfts_,$(xfts_bench_qrm_SOURCES))
dfts_bench_qrm_SOURCES := $(subst xfts_,dfts_,$(xfts_bench_qrm_SOURCES))
sfts_bench_qrm_SOURCES := $(subst xfts_,sfts_,$(xfts_bench_qrm_SOURCES))
zfts_bench_qrm_SOURCES := $(subst xfts_,zfts_,$(xfts_bench_qrm_SOURCES))

cfts_test_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(cfts_test_SOURCES)))
dfts_test_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(dfts_test_SOURCES)))
sfts_test_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(sfts_test_SOURCES)))
zfts_test_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(zfts_test_SOURCES)))

cfts_bench_qrm_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(cfts_bench_qrm_SOURCES)))
dfts_bench_qrm_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(dfts_bench_qrm_SOURCES)))
sfts_bench_qrm_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(sfts_bench_qrm_SOURCES)))
zfts_bench_qrm_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(zfts_bench_qrm_SOURCES)))


cfts_test_LDADD = $(subst xfts_,cfts_,$(xfts_test_LDADD))
dfts_test_LDADD = $(subst xfts_,dfts_,$(xfts_test_LDADD))
sfts_test_LDADD = $(subst xfts_,sfts_,$(xfts_test_LDADD))
zfts_test_LDADD = $(subst xfts_,zfts_,$(xfts_test_LDADD))

cfts_bench_qrm_LDADD = $(subst xfts_,cfts_,$(xfts_bench_qrm_LDADD))
dfts_bench_qrm_LDADD = $(subst xfts_,dfts_,$(xfts_bench_qrm_LDADD))
sfts_bench_qrm_LDADD = $(subst xfts_,sfts_,$(xfts_bench_qrm_LDADD))
zfts_bench_qrm_LDADD = $(subst xfts_,zfts_,$(xfts_bench_qrm_LDADD))


libcfts_src_a_SOURCES = $(subst xfts_,cfts_,$(libxfts_src_a_SOURCES))
libdfts_src_a_SOURCES = $(subst xfts_,dfts_,$(libxfts_src_a_SOURCES))
libsfts_src_a_SOURCES = $(subst xfts_,sfts_,$(libxfts_src_a_SOURCES))
libzfts_src_a_SOURCES = $(subst xfts_,zfts_,$(libxfts_src_a_SOURCES))

libcfts_src_a_OBJECTS = $(patsubst %.F90,%.o,$(filter %.F90,$(libcfts_src_a_SOURCES)))
libdfts_src_a_OBJECTS = $(patsubst %.F90,%.o,$(filter %.F90,$(libdfts_src_a_SOURCES)))
libsfts_src_a_OBJECTS = $(patsubst %.F90,%.o,$(filter %.F90,$(libsfts_src_a_SOURCES)))
libzfts_src_a_OBJECTS = $(patsubst %.F90,%.o,$(filter %.F90,$(libzfts_src_a_SOURCES)))


$(libcfts_src_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(QRM_FCFLAGS)
$(libdfts_src_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(QRM_FCFLAGS)
$(libsfts_src_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(QRM_FCFLAGS)
$(libzfts_src_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(QRM_FCFLAGS)

ifneq (,$(filter c,$(ARITHS)))
include $(cfts_test_OBJECTS:.o=.d)
endif

ifneq (,$(filter s,$(ARITHS)))
include $(sfts_test_OBJECTS:.o=.d)
endif

ifneq (,$(filter d,$(ARITHS)))
include $(dfts_test_OBJECTS:.o=.d)
endif

ifneq (,$(filter z,$(ARITHS)))
include $(zfts_test_OBJECTS:.o=.d)
endif


ifneq (,$(filter c,$(ARITHS)))
include $(cfts_bench_qrm_OBJECTS:.o=.d)
endif

ifneq (,$(filter s,$(ARITHS)))
include $(sfts_bench_qrm_OBJECTS:.o=.d)
endif

ifneq (,$(filter d,$(ARITHS)))
include $(dfts_bench_qrm_OBJECTS:.o=.d)
endif

ifneq (,$(filter z,$(ARITHS)))
include $(zfts_bench_qrm_OBJECTS:.o=.d)
endif


ifneq (,$(filter c,$(ARITHS)))
include $(libcfts_src_a_OBJECTS:.o=.d)
endif
#    
ifneq (,$(filter d,$(ARITHS)))
include $(libdfts_src_a_OBJECTS:.o=.d)
endif
#                                                                                                                                                                                                           
ifneq (,$(filter s,$(ARITHS)))
include $(libsfts_src_a_OBJECTS:.o=.d)
endif
#                                                                                                                                                                                                           
ifneq (,$(filter z,$(ARITHS)))
include $(libzfts_src_a_OBJECTS:.o=.d)
endif
