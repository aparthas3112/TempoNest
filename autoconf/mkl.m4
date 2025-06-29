dnl @synopsis SWIN_LIB_MKL
dnl 
AC_DEFUN([SWIN_LIB_MKL],
[
  AC_PROVIDE([SWIN_LIB_MKL])

  SWIN_PACKAGE_OPTIONS([mkl])

  AC_MSG_CHECKING([for Intel Math Kernel Library (MKL) installation])


  if test "$have_mkl" != "user disabled"; then

    SWIN_PACKAGE_FIND([mkl],[mkl_lapacke.h])
    SWIN_PACKAGE_TRY_COMPILE([mkl],[#include <mkl_lapacke.h>])

    SWIN_PACKAGE_TRY_LINK([mkl],[#include <mkl.h>],
                          [LAPACKE_dlamch('E');],
                          [-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm  -lmkl_avx2  -lmkl_def -qopenmp -liomp5])

  fi

  AC_MSG_RESULT([$have_mkl])


  if test x"$have_mkl" = xyes; then
      have_mkl=old
      AC_DEFINE([HAVE_MKL],[1],
                [Define if the old Intel Math Kernel Library is present])
      [$1]


      MKL_LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm  -lmkl_avx2  -lmkl_def -qopenmp -liomp5"
      MKL_CFLAGS="-DMKL_ILP64"
  else
      MKL_LIBS=""
      MKL_CFLAGS=""

  fi

  AC_SUBST(MKL_LIBS)
  AC_SUBST(MKL_CFLAGS)
  AM_CONDITIONAL(HAVE_MKL,[test "$have_mkl" = yes])

])

