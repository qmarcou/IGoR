ACLOCAL_AMFLAGS = -I ../m4

# Declare program name and install it in $(bindir) directory 
bin_PROGRAMS = igor

# List all Igor sources
igor_SOURCES = Aligner.cpp Aligner.h Bestscenarioscounter.cpp Bestscenarioscounter.h Counter.cpp Counter.h Coverageerrcounter.cpp Coverageerrcounter.h Deletion.cpp Deletion.h Dinuclmarkov.cpp Dinuclmarkov.h Errorrate.cpp Errorrate.h Genechoice.cpp Genechoice.h GenModel.cpp GenModel.h Hypermutationglobalerrorrate.cpp Hypermutationglobalerrorrate.h Insertion.cpp Insertion.h IntStr.cpp IntStr.h main.cpp Model_marginals.cpp Model_marginals.h Model_Parms.cpp Model_Parms.h Pgencounter.cpp Pgencounter.h Rec_Event.cpp Rec_Event.h Singleerrorrate.cpp Singleerrorrate.h Utils.cpp Utils.h

# Include GSL subparts and jemalloc without installation

igor_LDADD = ../libs/gsl_sub/libgsl_sub.la ../libs/jemalloc/lib/libjemalloc.a -lpthread #Need to link pthread since it is not statically linked with jemalloc (for now)

igor_CXXFLAGS = -I$(srcdir)/../libs/jemalloc/include/jemalloc $(OPENMP_CXXFLAGS)

#igor_AM_CPPFLAGS = -std=c++11 -I$(srcdir)/../libs/gsl_sub 
#igor_AM_CXXFLAGS = -I$(srcdir)/../libs/jemalloc/include/jemalloc $(OPENMP_CXXFLAGS)
#igor_AM_CFLAGS = $(OPENMP_CXXFLAGS)
# Compile using c++11 and include headers from $(srcdir)/../libs/gsl_sub to mimick the use of normal GSL and include it using #<gsl/gsl_foo.h>
#-I$(srcdir)/../libs/jemalloc doesn't really need to have the headers #-fopenMP would not work on every compiler, using AC_OPENMP macro instead
