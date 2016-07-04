# routines from the slatec library used in MaPHyS
libslatec_a_SOURCES = $(addprefix $(topsrcdir)/slatec/, \
	fdump.f i1mach.f icopy.f ipsort.f isort.f j4save.f \
	xercnt.f xerhlt.f xermsg.f xerprn.f xersve.f xgetua.f )
libslatec_a_OBJECTS = $(libslatec_a_SOURCES:.f=.o)
