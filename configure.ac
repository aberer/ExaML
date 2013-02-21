#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ExaML],[1.0.1],[bla@bla.org])
AC_CONFIG_SRCDIR([examl/axml.c])
AC_CONFIG_HEADERS([examl/config.h])

AM_INIT_AUTOMAKE(  examl  , 1.0.1 )

AM_SILENT_RULES([yes])



# Checks for programs.
AC_PROG_AWK
AC_PROG_CC
AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([pthread], [pthread_attr_init])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdint.h stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([getcwd gettimeofday memset pow sqrt strcasecmp strchr])

AX_EXT
AX_PROG_CC_MPI

AC_ARG_ENABLE(debug,
  [ --enable-debug    compile with debug information and  less optimization ],
  [case "${enableval}" in
     yes | no ) WITH_DEBUG="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
   esac],
  [WITH_DEBUG="no"]
)

AM_CONDITIONAL([WITH_DEBUG], [test  "x$WITH_DEBUG" = "xyes"])

AM_CONDITIONAL([WITH_AVX], [test "x$ax_cv_have_avx_ext" = "xyes" ] ) 

AC_ARG_ENABLE(productive,
  [ --disable-productive    developmental version of the code ],
  [case "${enableval}" in
     yes | no ) PRODUCTIVE="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-productive) ;;
   esac],
  [PRODUCTIVE="yes"]
)
AM_CONDITIONAL([PRODUCTIVE], [test  "x$PRODUCTIVE_DEBUG" = "xyes"])

AC_DEFINE([PRODUCTIVE], [0], ["is 1, if in productive mode "])
AS_IF([test "x$PRODUCTIVE" = "xyes"], AC_DEFINE(PRODUCTIVE, 1) )
AM_CONDITIONAL([PRODUCTIVE], [test  "x$PRODUCTIVE" = "xyes"])


AC_CONFIG_FILES([Makefile]) # TODO parser/Makefile 
AC_OUTPUT
