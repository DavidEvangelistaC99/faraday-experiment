#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
#  Main Madrigal configure.ac file
#
#  $Id: configure.ac 7143 2020-07-29 14:17:16Z brideout $

AC_PREREQ(2.59)
AC_INIT([Madrigal], [3.1.3], [brideout@haystack.mit.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_PROG_AR
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_ENABLE_SHARED
LT_INIT

# Checks for programs.
AC_PROG_CC
AC_PROG_F77

# Mix Fortran and C
AC_F77_LIBRARY_LDFLAGS
AC_F77_WRAPPERS

# Checks for libraries
 
# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF([long])
AM_CONDITIONAL([IS_32BIT], [test $ac_cv_sizeof_long = 4])

# Checks for library functions.

AC_CONFIG_FILES([
    Makefile
    madf/geolib/Makefile
    madf/applications/Makefile
    madc/madrec/Makefile
])
AC_OUTPUT
