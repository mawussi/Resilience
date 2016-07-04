# CERFACS libpackGMRES library shipped within MaPHyS
libpackgmres_a_SOURCES = $(addprefix $(topsrcdir)/packgmres/, \
	cPackgmres.f dPackgmres.f sPackgmres.f zPackgmres.f )
libpackgmres_a_OBJECTS = $(libpackgmres_a_SOURCES:.f=.o)
