ACLOCAL_AMFLAGS = -I ../m4

# Declare program name and install it in $(bindir) directory 
bin_PROGRAMS = igor 

# List all Igor sources
SOURCES = Aligner.cpp Aligner.h Bestscenarioscounter.cpp Bestscenarioscounter.h CDR3SeqData.h CDR3SeqData.cpp Counter.cpp Counter.h Coverageerrcounter.cpp Coverageerrcounter.h Deletion.cpp Deletion.h Dinuclmarkov.cpp Dinuclmarkov.h Errorscounter.cpp Errorscounter.h Errorrate.cpp Errorrate.h ExtractFeatures.h ExtractFeatures.cpp Genechoice.cpp Genechoice.h GenModel.cpp GenModel.h HypermutationfullNmererrorrate.cpp HypermutationfullNmererrorrate.h Hypermutationglobalerrorrate.cpp Hypermutationglobalerrorrate.h Insertion.cpp Insertion.h IntStr.cpp IntStr.h Model_marginals.cpp Model_marginals.h Model_Parms.cpp Model_Parms.h Pgencounter.cpp Pgencounter.h Rec_Event.cpp Rec_Event.h Singleerrorrate.cpp Singleerrorrate.h Utils.cpp Utils.h

igor_SOURCES = $(SOURCES) main.cpp


# Include GSL subparts and jemalloc without installation

igor_LDADD = ../libs/gsl_sub/libgsl_sub.la ../libs/jemalloc/lib/libjemalloc.a -lpthread -ldl #Need to link pthread since it is not statically linked with jemalloc (for now)


igor_CXXFLAGS = -I$(srcdir)/../libs/jemalloc/include/jemalloc -I$(srcdir)/../libs/gsl_sub $(OPENMP_CXXFLAGS) -DIGOR_DATA_DIR=\"$(pkgdatadir)\"


# Pass IGoR data install directory path through the compiler with -DIGOR_DATA_DIR=\"$(pkgdatadir)\"
# and not in config.h (as it is best practice according to: https://stackoverflow.com/questions/5867136/autoconf-how-to-get-installation-paths-into-config-h#11317815)



#igor_AM_CPPFLAGS = -std=c++11 -I$(srcdir)/../libs/gsl_sub 
#igor_AM_CXXFLAGS = -I$(srcdir)/../libs/jemalloc/include/jemalloc $(OPENMP_CXXFLAGS)
#igor_AM_CFLAGS = $(OPENMP_CXXFLAGS)
# Compile using c++11 and include headers from $(srcdir)/../libs/gsl_sub to mimick the use of normal GSL and include it using #<gsl/gsl_foo.h>
#-I$(srcdir)/../libs/jemalloc doesn't really need to have the headers #-fopenMP would not work on every compiler, using AC_OPENMP macro instead

