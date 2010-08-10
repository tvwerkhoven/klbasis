## @file configure.ac
## @author Tim van Werkhoven (t.i.m.vanwerkhoven@xs4all.nl)

AC_INIT([klbasis], [0.1], [t.i.m.vanwerkhoven@xs4all.nl])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_FILES([
Makefile
])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_RANLIB

# These libraries are crucial for all targets
AC_SEARCH_LIBS([sin],
		[m],
		[],
		[AC_MSG_ERROR([Cannot build libraries without libm])])

AC_SEARCH_LIBS([pthread_create],
		[pthread],
		[],
		[AC_MSG_ERROR([Cannot build libraries without pthreads])])

PKG_CHECK_MODULES(GSL, [gsl], [],
		[AC_MSG_ERROR([Cannot build libraries without gsl])])

AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_OUTPUT
