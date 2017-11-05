#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([igor], [1.0.0.2], [quentin.marcou@lpt.ens.fr])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

#AC_CONFIG_SRCDIR([igor_src/Deletion.h])
AC_CONFIG_HEADERS([config.h])

#INIT LIBTOOL
AC_CONFIG_MACRO_DIRS([m4])
LT_INIT

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

#Force usage of C++ and OPENMP
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_11
AC_OPENMP
#AC_SUBST(AM_CXXFLAGS,$OP‌​ENMP_CXXFLAGS)

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 igor_src/Makefile
                 libs/Makefile])
AC_CONFIG_SUBDIRS([libs/gsl_sub
                   libs/jemalloc
					igor_src])
AC_OUTPUT
