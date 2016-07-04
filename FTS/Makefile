SHELL=/bin/bash
topsrcdir   := $(shell pwd)
abstopsrcdir:= $(abspath $(topsrcdir))
buildmoddir := $(strip $(abstopsrcdir)/modules)
buildlibdir := $(strip $(topsrcdir)/lib )
scriptdir   := $(topsrcdir)/scripts

include ${topsrcdir}/Makefile.inc

##
##  General inference rules.
##




ARITHS ?=c d s z


# Target 

lib_LIBRARIES:=$(addprefix $(buildlibdir)/, \
	libftsolver.a libpackgmres.a libpackfgmres.a  libslatec.a)

test_LIBRARIES += $(addprefix $(buildlibdir)/, $(ARITHS:%=lib%fts_toolkit.a))
bench_qrm_LIBRARIES += $(addprefix $(buildlibdir)/, $(ARITHS:%=lib%fts_toolkit.a))


test_PROGRAMS := $(addprefix $(topsrcdir)/src/, $(ARITHS:%=%fts_test) )
bench_qrm_PROGRAMS := $(addprefix $(topsrcdir)/src/, $(ARITHS:%=%fts_bench_qrm) )



.SECONDEXPANSION:
.PHONY : all

.PHONY : clean
.PHONY : test
.PHONY : bench_qrm

all : $(lib_LIBRARIES) 

test : $(test_PROGRAMS)
bench_qrm : $(bench_qrm_PROGRAMS)



$(test_PROGRAMS) : $(topsrcdir)/src/%$(EXEEXT) : $$(%_OBJECTS) $$(%_LDADD)
	$(FC) -o $@ $(filter %.o %.a,$^) $(ALL_LDFLAGS) $(EXTRA_LDFLAGS) 

$(test_PROGRAMS) : $(test_LIBRARIES)

$(bench_qrm_PROGRAMS) : $(topsrcdir)/src/%$(EXEEXT) : $$(%_OBJECTS) $$(%_LDADD)
	$(FC) -o $@ $(filter %.o %.a,$^) $(ALL_LDFLAGS) $(EXTRA_LDFLAGS) 

$(bench_qrm_PROGRAMS) : $(bench_qrm_LIBRARIES)


## Build a object
%$(OBJEXT) : %.F90
	cd $(buildmoddir) ; $(FC)  $(ALL_FCFLAGS) $(EXTRA_FCFLAGS) -o $(abspath $@) -c $(abspath $<)

%$(OBJEXT) : %.F
	cd $(buildmoddir) ; $(F77)  $(ALL_FCFLAGS) $(EXTRA_FFLAGS) -o $(abspath $@) -c $(abspath $<)	
#	$(F77) $(ALL_FCFLAGS) $(EXTRA_FFLAGS) -o $@ -c $<

%$(OBJEXT) : %.f
	$(F77)  $(ALL_FFLAGS) $(EXTRA_FFLAGS) -o $@ -c $<

%$(OBJEXT) : %.c
	$(CC)  $(ALL_CFLAGS) $(EXTRA_CFLAGS) -o $@ -c $<


## Build a static library
$(buildlibdir)/%.a : $$(%_a_OBJECTS) $$(%_a_LIBADD)
	$(AR) $(ARFLAGS) $@ $?

## Build an arithmetic dependent source file

x2c=$(scriptdir)/arithmetic_preprocessing c
x2d=$(scriptdir)/arithmetic_preprocessing d
x2s=$(scriptdir)/arithmetic_preprocessing s
x2z=$(scriptdir)/arithmetic_preprocessing z

cfts_%.F90 : xfts_%.F90;	$(x2c) $< $@ 
dfts_%.F90 : xfts_%.F90;	$(x2d) $< $@ 
sfts_%.F90 : xfts_%.F90;	$(x2s) $< $@ 
zfts_%.F90 : xfts_%.F90;	$(x2z) $< $@ 

cfts_%.F90 : xfts_%.F;	$(x2c) $< $@ 
dfts_%.F90 : xfts_%.F;	$(x2d) $< $@ 
sfts_%.F90 : xfts_%.F;	$(x2s) $< $@ 
zfts_%.F90 : xfts_%.F;	$(x2z) $< $@ 



### Build the dependency files

# copy .mod of qr_mumps in modules


ifeq (,$(findstring YES,$(WEAKLY_LINK_OBJECTS2MODULES))) 
DEPENDS_F90_OPTS := -p $(buildmoddir)/
else
DEPENDS_F90_OPTS := -w -p $(buildmoddir)/ 
endif

%.d : %.F90
	$(scriptdir)/depends_f90.sh $(DEPENDS_F90_OPTS) $< > $@

%.d : %.F
	$(scriptdir)/depends_f90.sh -f ".F" $(DEPENDS_F90_OPTS) $< > $@

%.d : %.f
	touch $@

%.d : %.c
	touch $@



ifeq (,$(filter %clean,$(MAKECMDGOALS))) # avoid unecessary inclusion  
include ${topsrcdir}/slatec/rules.mk
include ${topsrcdir}/packgmres/rules.mk
include ${topsrcdir}/packfgmres/rules.mk
include ${topsrcdir}/part/rules.mk
include ${topsrcdir}/sm/rules.mk
include ${topsrcdir}/dm/rules.mk
include ${topsrcdir}/common/rules.mk
include ${topsrcdir}/lib/rules.mk
include ${topsrcdir}/toolkit/rules.mk
include ${topsrcdir}/src/rules.mk
endif


## Build a module
%.o:%.mod # cancel Modula-2 rules which conflict with Fortran rules
$(buildmoddir)/%$(MODEXT):$(topsrcdir)/src/%$(OBJEXT) ;
$(buildmoddir)/%$(MODEXT):$(topsrcdir)/common/%$(OBJEXT) ;
$(buildmoddir)/%$(MODEXT):$(topsrcdir)/dm/%$(OBJEXT) ;	   
$(buildmoddir)/%$(MODEXT):$(topsrcdir)/sm/%$(OBJEXT) ;	  
$(buildmoddir)/%$(MODEXT):$(topsrcdir)/part/%$(OBJEXT) ;  
$(buildmoddir)/%$(MODEXT):$(topsrcdir)/toolkit/%$(OBJEXT) ;



.PRECIOUS : %$(MODEXT)  # generated modules
.PRECIOUS : %.d
.PRECIOUS : %.F90

ifneq (,$(findstring YES,$(KEEP_GENERATED_SOURCES))) # DO NOT WORK WHY ?


.SECONDARY : $(libfts_common_a_SOURCES) 
.SECONDARY : \
	$(libcfts_dm_a_SOURCES) $(libcfts_sm_a_SOURCES) \
	$(libcfts_toolkit_a_SOURCES) $(libcfts_part_a_SOURCES) \
	$(cfts_test_SOURCES) $(cfts_bench_qrm_SOURCES) $(libcfts_src_a_SOURCES) 

.SECONDARY : \
	$(libdfts_dm_a_SOURCES) $(libdfts_sm_a_SOURCES) \
	$(libdfts_toolkit_a_SOURCES) $(libdfts_part_a_SOURCES) \
	$(dfts_test_SOURCES)  $(dfts_bench_qrm_SOURCES) $(libdfts_src_a_SOURCES) 


.SECONDARY : \
	$(libsfts_dm_a_SOURCES) $(libsfts_sm_a_SOURCES) \
	$(libsfts_toolkit_a_SOURCES) 	$(libsfts_part_a_SOURCES) \
	$(sfts_test_SOURCES) $(sfts_bench_qrm_SOURCES) $(libsfts_src_a_SOURCES) 


.SECONDARY : \
	$(libzfts_dm_a_SOURCES) $(libzfts_sm_a_SOURCES) \
	$(libzfts_toolkit_a_SOURCES)	$(libzfts_part_a_SOURCES) \
	$(zfts_test_SOURCES) $(zfts_bench_qrm_SOURCES) $(libzfts_src_a_SOURCES) 

endif



SUBDIRS:=common dm sm   part src \
         toolkit slatec  packgmres packfgmres


## Generate the tags
.PHONY: TAGS
TAGS : 
	-etags $(foreach d,$(SUBDIRS) include,\
	$(wildcard $(d)/*.f $(d)/*.F $(d)/*.c $(d)/*.h \
	$(d)/xfts_*.F90 $(d)/fts_*.F90 ))


clean : mostlyclean
	$(foreach d,$(SUBDIRS),cd $(abstopsrcdir)/$(d); $(RM) *.d [cdsz]fts_*.F{90,} ;)

mostlyclean :
	$(RM) $(buildlibdir)/*.a *~
	cd $(abstopsrcdir)/src; $(RM) [cdsz]fts_test [cdsz]fts_bench_qrm
	cd $(abstopsrcdir)/modules; $(RM) *
	$(foreach d,$(SUBDIRS),cd $(abstopsrcdir)/$(d); $(RM) *$(OBJEXT) *$(MODEXT) *~ ;)
	$(RM) $(buildmoddir)/*$(MODEXT)
#	cp ${QRM_prefix}/src/modules/*mod  modules/			
##
##  Project rules.
##
