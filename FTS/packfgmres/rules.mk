# CERFACS libpackFGMRES library shipped within MaPHyS
libpackfgmres_a_SOURCES = $(addprefix $(topsrcdir)/packfgmres/, \
	cPackfgmres.f dPackfgmres.f sPackfgmres.f zPackfgmres.f )
libpackfgmres_a_OBJECTS = $(libpackfgmres_a_SOURCES:.f=.o)
