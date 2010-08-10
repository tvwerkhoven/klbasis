AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -g -Wall -O3 -ftree-vectorize 
AM_CPPFLAGS = $(GSL_CFLAGS) -DHAVE_INLINE \
		-DGSL_RANGE_CHECK_OFF
AM_LDFLAGS = $(GSL_LIBS)

# Program directives
bin_PROGRAMS = klbasis
klbasis_SOURCES = a1a2.cc kernel.cc klbasis.cc util.cc
klbasis_SOURCES += a1a2.h kernel.h klbasis.h util.h
