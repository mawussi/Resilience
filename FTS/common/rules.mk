# maphys common library (common means common to all arithmetics)
libfts_common_a_SOURCES = $(addprefix $(topsrcdir)/common/,\
	fts_common.F90 \
	fts_dbg_mod.F90 fts_env_mod.F90 \
	fts_error_mod.F90 fts_log_mod.F90 fts_mem_mod.F90 \
	fts_time_mod.F90)
libfts_common_a_OBJECTS = $(libfts_common_a_SOURCES:%.F90=%$(OBJEXT))

### 
include $(libfts_common_a_OBJECTS:.o=.d)