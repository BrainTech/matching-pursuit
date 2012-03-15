#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([MP Engine], [5],
        [rkus@fuw.edu.pl],
        [mp_engine], [https://escher.fuw.edu.pl/git/mp_engine])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE(foreign dist-bzip2 dist-xz std-options silent-rules
		 subdir-objects tar-pax -Wall -Werror -Wno-portability)
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_SRCDIR([src/mp5.c])
AC_CONFIG_HEADERS([config/config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# This makes sure pkg.m4 is available.
m4_pattern_forbid([^_?PKG_[A-Z_]+$],
	[*** pkg.m4 missing, please install pkg-config])

AC_SEARCH_LIBS([atan2], [m])
AC_SEARCH_LIBS([sincos], [m])

AC_ARG_WITH([fftw3],
	AS_HELP_STRING([--with-fftw3=...], [explicit path to FFTW3 DLL]))
AS_CASE([$with_fftw3],
	[yes],
	[AC_MSG_ERROR([--with-fftw3 requires an explicit path])],
	[no],
	[AC_MSG_ERROR([--without-fftw3 is meaningless])],
	[''],
	[],
	[[FFTW3_LIBS=$with_fftw3]
	 [FFTW3_CFLAGS=-I`dirname $with_fftw3`]
	 AC_MSG_RESULT([FFTW3 library detection overriden])
	 AC_SUBST([FFTW3DLL], [$with_fftw3])
	]
)

PKG_CHECK_MODULES(FFTW3, [ fftw3 ])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
