## Process this file with automake to produce Makefile.in

# AUTOMAKE_OPTIONS = readme-alpha

SUBDIRS = gsl utils sys test err block vector matrix permutation cblas blas linalg rng histogram

SUBLIBS = block/libgslblock.la blas/libgslblas.la cblas/libgslcblas.la err/libgslerr.la histogram/libgslhistogram.la linalg/libgsllinalg.la matrix/libgslmatrix.la permutation/libgslpermutation.la rng/libgslrng.la test/libgsltest.la utils/libutils.la vector/libgslvector.la

ACLOCAL_AMFLAGS = -I ../../m4

# pkginclude_HEADERS = gsl_math.h gsl_pow_int.h gsl_nan.h gsl_machine.h gsl_mode.h gsl_precision.h gsl_types.h gsl_version.h gsl_minmax.h gsl_inline.h
#This line was commented out since we do not wish to install these GSL subparts

noinst_SCRIPTS = gsl-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA= gsl.pc

CLEANFILES = gsl.pc gsl-config

EXTRA_DIST = autogen.sh gsl-config.in gsl.pc.in configure.ac THANKS BUGS gsl.spec.in gsl.m4 test_gsl_histogram.sh pkgconfig.test

#noinst_LIBRARIES = libgsl_sub.a #Do not install, does not work with .la???
noinst_LTLIBRARIES = libgsl_sub.la
libgsl_sub_la_SOURCES = version.c
libgsl_sub_la_LIBADD = $(GSL_LIBADD) $(SUBLIBS)
libgsl_sub_la_LDFLAGS = $(GSL_LDFLAGS) -version-info $(GSL_LT_VERSION)
noinst_HEADERS = templates_on.h templates_off.h build.h

m4datadir = $(datadir)/aclocal
m4data_DATA = gsl.m4

#noinst_PROGRAMS = gsl-randist gsl-histogram
noinst_PROGRAMS = gsl-histogram #Do not install but still build in order to produce the tests

#gsl_randist_SOURCES = gsl-randist.c
#gsl_randist_LDADD = libgsl_sub.la cblas/libgslcblas.la

gsl_histogram_SOURCES = gsl-histogram.c
gsl_histogram_LDADD = libgsl_sub.la #cblas/libgslcblas.la compiled straight with the rest of the lib

check_SCRIPTS = test_gsl_histogram.sh pkgconfig.test
TESTS = test_gsl_histogram.sh pkgconfig.test

#bin_PROGRAMS = main dummy
#dummy_SOURCES = version.c
#dummy_LDADD = $(SUBLIBS)
#main_SOURCES = version.c env.c
#main_LDADD = libgsl.la

edit = $(SED) \
	-e 's|@prefix[@]|$(prefix)|g' \
	-e 's|@exec_prefix[@]|$(exec_prefix)|g' \
	-e 's|@libdir[@]|$(libdir)|g' \
	-e 's|@includedir[@]|$(includedir)|g' \
	-e 's|@GSL_CFLAGS[@]|$(GSL_CFLAGS)|g' \
	-e 's|@GSL_LIBM[@]|$(GSL_LIBM)|g' \
	-e 's|@GSL_LIBS[@]|$(GSL_LIBS)|g' \
	-e 's|@LIBS[@]|$(LIBS)|g' \
	-e 's|@VERSION[@]|$(VERSION)|g'

gsl-config gsl.pc: Makefile 
	@rm -f $@ $@.tmp
	@$(edit) '$(srcdir)/$@.in' >>$@.tmp
	@chmod a-w $@.tmp
	@mv $@.tmp $@
	@echo creating $@

gsl-config: $(srcdir)/gsl-config.in
gsl.pc: $(srcdir)/gsl.pc.in

