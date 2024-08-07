#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([TempoNest],[0.1],[ltl21@mrao.cam.ac.uk])

# Check for MPI
AC_LANG_PUSH([C++])
AC_LANG_POP([C++])

# Set MPI flags
AC_SUBST(MPI_CXXFLAGS)
AC_SUBST(MPI_CXXLIBS)
AC_SUBST(MPI_CXXLDFLAGS)

# Capture MPI compiler flags
if test "x$MPICXX" != "x"; then
    MPI_CXXFLAGS=`$MPICXX -showme:compile`
    MPI_CXXLIBS=`$MPICXX -showme:link`
    MPI_CXXLDFLAGS=`$MPICXX -showme:link | sed 's/-l[^ ]*//g'`
fi

AC_ARG_VAR([TEMPO2], [Default TEMPO2 installation directory])
AC_PREFIX_DEFAULT([${TEMPO2:-"/usr/local"}])

if test -n "${CXXFLAGS}"; then
	user_set_cxxflags=yes
fi

if test -n "${CFLAGS}"; then
	user_set_cflags=yes
fi

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([autoconf])

# Checks for programs.
AC_PROG_CXX([mpicxx])
AC_PROG_CXXCPP([mpicxx])
AC_PROG_CC([mpicc])
AC_PROG_CPP([mpicxx])
AM_PROG_CC_C_O([mpicc])
AC_PROG_F77([gfortran])
AC_PROG_F77_C_O([gfortran])
AC_F77_LIBRARY_LDFLAGS
AC_PROG_FC

AC_PROG_INSTALL
AC_PROG_LIBTOOL

# Optimization flags
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [Enable debug build])],
    [enable_debug=$enableval],
    [enable_debug=no])

AS_IF([test "x$enable_debug" = xyes], [
    CFLAGS="-g -O0"
    CXXFLAGS="-g -O0"
    FFLAGS="-g -O0 -DHAVE_CONFIG_H"
], [
    CFLAGS="-O3 -march=native -mtune=native -ffast-math"
    CXXFLAGS="-O3 -march=native -mtune=native -ffast-math"
    FFLAGS="-O3 -march=native -mtune=native -ffast-math -DHAVE_CONFIG_H"
])

# Add warning flags
CFLAGS="$CFLAGS -Wall -Wextra -Wpedantic"
CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wpedantic"

# Add MPI flags
CXXFLAGS="$CXXFLAGS $MPI_CXXFLAGS"
LDFLAGS="$LDFLAGS $MPI_CXXLDFLAGS"

# Respect user-set flags
if test "x$user_set_cxxflags" = xyes; then
    AC_MSG_NOTICE([Using user-provided CXXFLAGS: $CXXFLAGS])
fi
if test "x$user_set_cflags" = xyes; then
    AC_MSG_NOTICE([Using user-provided CFLAGS: $CFLAGS])
fi
if test "x$user_set_fflags" = xyes; then
    AC_MSG_NOTICE([Using user-provided FFLAGS: $FFLAGS])
fi

AC_CONFIG_AUX_DIR([config])
#AC_CONFIG_SRCDIR([Makefile.am])

# Enable large-file support
AC_SYS_LARGEFILE

# Determine the machine endian
AC_C_BIGENDIAN

if test $can_build_shared == no ; then
	AC_WARN("Libtool cannot build shared libraries - ")
fi

# Check whether the 'cos' function can be found in library 'm'
# (standard C math library).
# Math library is required for the tempo2 predictor library
AC_CHECK_LIB(m, cos)

# look for third-party packages in various places (or not)
#SWIN_LOCAL

AC_MSG_CHECKING([for \$TEMPO2])
if test -n "$TEMPO2" ; then
	TEMPO2_DIR=$TEMPO2
	AC_SUBST(TEMPO2_DIR)
	have_tempo2_dir=yes
else
	have_tempo2_dir=no
fi

AC_MSG_RESULT($have_tempo2_dir)

if test X$have_tempo2_dir != Xyes; then
	AC_WARN("You must set the $TEMPO2 variable before configuring.")
fi

AC_ARG_WITH([tempo2-plug-dir],
              AC_HELP_STRING([--with-tempo2-plug-dir=DIR],
                             [Specify DIR to install plugins]))


AC_MSG_CHECKING([for plugin install path])
if test -z "$with_tempo2_plug_dir" ; then
	TEMPO2_PLUG_INSTALL=$TEMPO2/plugins/
else
	TEMPO2_PLUG_INSTALL=$with_tempo2_plug_dir
	AC_DEFINE_UNQUOTED([TEMPO2_CONFIGURE_PLUG], "$with_tempo2_plug_dir",["Tempo2 plugin path"])
fi

AC_SUBST(TEMPO2_PLUG_INSTALL)
AC_MSG_RESULT($TEMPO2_PLUG_INSTALL)

MJK_FORTRAN_TEST

if test X$fortran_c_links != Xyes ; then
	AC_MSG_ERROR(Your C and F77 compilers do not cross link (\$CC=$CC \$F77=$F77))
fi

LT_LIB_DLLOAD

#SWIN_LIB_PGPLOT
#SWIN_LIB_FFTW
#SWIN_LIB_CFITSIO



#GSL basic stuff

AX_PATH_GSL(1.11.0, have_gsl=yes,have_gsl=no)
AM_CONDITIONAL(HAVE_GSL, [test x"$have_gsl" = xyes])

AC_F77_WRAPPERS

ACX_BLAS(have_blas=yes,have_blas=no)
AC_MSG_CHECKING([for BLAS])
AC_MSG_RESULT($have_blas)
AM_CONDITIONAL(HAVE_BLAS, [test x$have_blas = xyes])

# Check for Eigen
AC_ARG_WITH([eigen],
    [AS_HELP_STRING([--with-eigen=DIR],
        [prefix for Eigen headers])],
    [EIGEN_CPPFLAGS="-I$withval"],
    [EIGEN_CPPFLAGS='-I/usr/include/eigen3'])

AC_MSG_CHECKING([for Eigen])
save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $EIGEN_CPPFLAGS"
AC_LANG_PUSH([C++])
AC_TRY_COMPILE([#include <Eigen/Dense>],
               [Eigen::MatrixXd m(2,2);],
               [have_eigen=yes],
               [have_eigen=no])
AC_LANG_POP([C++])
CPPFLAGS=$save_CPPFLAGS
AC_MSG_RESULT([$have_eigen])

if test "x$have_eigen" = xno; then
    AC_MSG_ERROR([
    Eigen is required but not found. Please install Eigen or specify its location:
    ./configure CPPFLAGS="-I/path/to/eigen"
    or
    ./configure --with-eigen=/path/to/eigen
    
    Current EIGEN_CPPFLAGS: $EIGEN_CPPFLAGS
    Current compile command: $CXX $CPPFLAGS $CXXFLAGS
    ])
fi

# Add Eigen CPPFLAGS to the global CPPFLAGS
CPPFLAGS="$CPPFLAGS $EIGEN_CPPFLAGS"
AC_SUBST([EIGEN_CPPFLAGS])

# Check for OpenMP
AC_OPENMP
if test "x$ac_cv_openmp" != "xno"; then
    OPENMP_CXXFLAGS="$OPENMP_CXXFLAGS"
else
    AC_MSG_WARN([OpenMP not found. Eigen will run single-threaded.])
fi
AC_SUBST(OPENMP_CXXFLAGS)

# test for LAPACK acceleration
ACX_LAPACK(have_lapack=yes, have_lapack=no)
AC_MSG_CHECKING([for LAPACK])
AC_MSG_RESULT($have_lapack)
AM_CONDITIONAL(HAVE_LAPACK, [test x"$have_lapack" = xyes])

SWIN_LIB_MULTINEST
SWIN_LIB_CULA
SWIN_LIB_TEMPO2
SWIN_LIB_MKL
#SWIN_LIB_PSRCHIVE
#SWIN_LIB_QDINSTALL
SWIN_LIB_MLAPACK

AC_MSG_CHECKING([for LOGIN_ARCH variable])
if test -n "$LOGIN_ARCH" ; then
	t2arch=$LOGIN_ARCH
	else if test -n "$OSTYPE" ; then
		AC_MSG_RESULT(no)
		AC_MSG_CHECKING([for OSTYPE variable])
		t2arch=$OSTYPE
		else if test -n $OSTYPE ; then
		AC_MSG_RESULT(no)
		AC_MSG_CHECKING([for uname])
			t2arch=`uname`
		fi
	fi
fi


AC_MSG_RESULT($t2arch)

AC_DEFINE_UNQUOTED([TEMPO2_ARCH], "$t2arch",["Tempo2 ARCH variable"])
T2ARCH=$t2arch

AC_SUBST(T2ARCH)


AM_INIT_AUTOMAKE([foreign subdir-objects])

AM_CONFIG_HEADER([config.h])
AC_CONFIG_FILES([Makefile])

#AC_CONFIG_FILES([Makefile plugin/Makefile sofa/Makefile
#		 unsupported_plugins/Makefile
#		 unsupported_plugins/toasim/Makefile
#		 unsupported_plugins/ppta/Makefile
#		 unsupported_plugins/glast/Makefile])

if test "x$have_eigen" = xyes; then
    AC_MSG_NOTICE([Eigen found and will be used])
else
    AC_MSG_NOTICE([Eigen not found. Some matrix operations may be slower])
fi

AC_MSG_NOTICE([Summary of build options:])
AC_MSG_NOTICE([  Debug mode: $enable_debug])
AC_MSG_NOTICE([  CFLAGS: $CFLAGS])
AC_MSG_NOTICE([  CXXFLAGS: $CXXFLAGS])
AC_MSG_NOTICE([  FFLAGS: $FFLAGS])
AC_MSG_NOTICE([  Eigen include flags: $EIGEN_CPPFLAGS])

AC_OUTPUT

AC_MSG_NOTICE(Fortran compiler is $F77 $FFLAGS)
AC_MSG_NOTICE(C compiler is $CC $CFLAGS)
AC_MSG_NOTICE(C++ compiler is $CXX $CXXFLAGS)
AC_MSG_NOTICE(Global linking flags are $LDFLAGS)


# User warnings

if test X$have_gsl != Xyes; then
	AC_MSG_NOTICE(Plugins: GSL was not found - plugins requiring GSL will not be compiled)
fi



if test X$have_lapack = Xyes; then
	AC_MSG_NOTICE(Using LAPACK library for Cholesky decomposition)
fi

if test X$have_blas = Xyes; then
	AC_MSG_NOTICE(Using BLAS library for matrix multiplication)
fi

if test X$have_multinest != Xyes; then
    if X$have_cula = Xyes; then
        AC_MSG_NOTICE(TempoNest will be built with CUDA/CULA GPU support)
    else
        AC_MSG_NOTICE(TempoNest will be built without GPU support)
    fi
fi

if test X$have_psrchive = Xyes; then
	AC_MSG_NOTICE(Using PSRCHIVE library)
else
	AC_MSG_NOTICE([PSRCHIVE library not found, will not compile ToA estimation part of TempoNest])
fi

#if test X$have_qdinstall = Xyes; then
#        AC_MSG_NOTICE(Using QD library)
#	if test X$have_mlapack = Xyes; then
#		AC_MSG_NOTICE(Using MLAPACK library)
#	else
#		AC_MSG_NOTICE(MLAPACK library not found, will not compile high precision part of TempoNest)
#	fi
#else
#        AC_MSG_NOTICE(QD library not found, will not compile high precision part of TempoNest)
#fi


AC_MSG_NOTICE(Tempo2 includes: $TEMPO2_CFLAGS)


if test X$can_build_shared = Xno ; then
	echo ""
	AC_WARN("Libtool cannot build shared libraries - no plugins can be built!")
	echo ""
else
	echo ""
	echo "NOTICE: TempoNest can now be built with 'make temponest && make temponest-install'"
fi


