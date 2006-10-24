
dnl Jasjeet S. Sekhon <sekhon@berkeley.edu> http://sekhon.berkeley.edu
dnl 2006-09-17
dnl Process this file with autoconf to produce a configure script.

dnl Matching Version 3.6-0+ we switched to configure.ac instead of the
dnl configure hack.

dnl Matching Version 3.4.4+: we need to add the -ffriend-injection
dnl argument for g++ versions 4.1+ because ARM-style name-injection of
dnl friend declarations is no longer the default.  We check for the
dnl flag directly instead of relying on gcc --version information because of 
dnl backporting issues---e.g., Gentoo.

dnl We also setup optimization flags.  This VASTLY improves performance.

dnl We use makevars on all platforms instead of
dnl the makefile because of build issues on OS X.  This requires moving
dnl files around because we do not want our cblas wrappers compiled on
dnl OS X, but we need them for all other platforms.

dnl On Windows configure does *not* need to be run.  On Windows,
dnl Makevars.win should simply be used because we can assume that we are
dnl using the mingw32 compiler

dnl -fstrict-aliasing -freorder-blocks -fsched-interblock is added for OS X 
dnl Already included in the Linux gcc -O2 option, but does no harm on Linux

AC_INIT(Matching, 4.0.8+, sekhon@berkeley.edu)

dnl MYCFLAGS0 should work on most compilers and MYCFLAGS1 should work on most GNU C compilers.  
dnl MYCFLAGS2 may not (generally needs gcc 4.0+).   The last three options in MYCFLAGS2 are the OS X 
dnl addons referenced above.



MYCFLAGS0="-O3"
MYCFLAGS1="-finline-functions -ffast-math -funroll-loops -fexpensive-optimizations"
MYCFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock"

dnl -funsafe-loop-optimizations for gcc 4.2+; the following doesn't work any better in any case
dnl MYCFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock -mfpmath=sse,387 -fgcse -fgcse-after-reload -fgcse-lm -fgcse-sm -fsched-spec-load-dangerous -funsafe-loop-optimizations"



MYCXXFLAGS0="-O3"
MYCXXFLAGS1="-finline-functions -ffast-math -funroll-loops -fexpensive-optimizations"
MYCXXFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock"

dnl -funsafe-loop-optimizations for gcc 4.2+; the following doesn't work any better in any case
dnl MYCXXFLAGS2="-funswitch-loops -fgcse-after-reload -fstrict-aliasing -freorder-blocks -fsched-interblock -mfpmath=sse,387 -fgcse -fgcse-after-reload -fgcse-lm -fgcse-sm -fsched-spec-load-dangerous -funsafe-loop-optimizations"


AC_PROG_CC

dnl Find the compilers to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`

AC_MSG_CHECKING([whether $CC accepts $MYCFLAGS0])
AC_LANG([C])
CFLAGS=$MYCFLAGS0
AC_COMPILE_IFELSE([[char b[10];]], [mycflags0_set=yes], [mycflags0_set=no])    
AC_MSG_RESULT($mycflags0_set)    

if test $mycflags0_set = no; then    
   PKG_CLFAGS=""
else
   PKG_CFLAGS=$MYCFLAGS0
fi

AC_MSG_CHECKING([whether $CC accepts $MYCFLAGS1])
AC_LANG([C])
CFLAGS=$MYCFLAGS1
AC_COMPILE_IFELSE([[char b[10];]], [mycflags1_set=yes], [mycflags1_set=no])    
AC_MSG_RESULT($mycflags1_set)    

if test $mycflags1_set = yes; then    
   PKG_CFLAGS="$PKG_CFLAGS $MYCFLAGS1"
fi

AC_MSG_CHECKING([whether $CC accepts $MYCFLAGS2])
AC_LANG([C])
CFLAGS=$MYCFLAGS2
AC_COMPILE_IFELSE([[char b[10];]], [mycflags2_set=yes], [mycflags2_set=no])    
AC_MSG_RESULT($mycflags2_set)    

if test $mycflags2_set = yes; then    
   PKG_CFLAGS="$PKG_CFLAGS $MYCFLAGS2"
fi

AC_PROG_CXX

if test $GXX = no; then
    AC_MSG_WARN([**** GNU g++ IS HIGHLY RECOMMENDED.  I will continue, but this may not work. ****])
fi

AC_MSG_CHECKING([whether $CXX accepts $MYCXXFLAGS0])
AC_LANG([C++])
CXXFLAGS=$MYCXXFLAGS0
AC_COMPILE_IFELSE([[char b[10];]], [mycxxflags0_set=yes], [mycxxflags0_set=no])    
AC_MSG_RESULT($mycxxflags0_set)    

if test $mycxxflags0_set = no; then    
  PKG_CXXLFAGS=""
else
  PKG_CXXFLAGS=$MYCXXFLAGS0
fi

AC_MSG_CHECKING([whether $CXX accepts $MYCXXFLAGS1])
AC_LANG([C++])
CXXFLAGS=$MYCXXFLAGS1
AC_COMPILE_IFELSE([[char b[10];]], [mycxxflags1_set=yes], [mycxxflags1_set=no])    
AC_MSG_RESULT($mycxxflags1_set)    

AC_MSG_CHECKING([whether $CXX accepts $MYCXXFLAGS2])
AC_LANG([C++])
CXXFLAGS=$MYCXXFLAGS2
AC_COMPILE_IFELSE([[char b[10];]], [mycxxflags2_set=yes], [mycxxflags2_set=no])    
AC_MSG_RESULT($mycxxflags2_set)    

if test $mycxxflags1_set = yes; then    
   PKG_CXXFLAGS="$PKG_CXXFLAGS $MYCXXFLAGS1"
fi

if test $mycxxflags2_set = yes; then    
   PKG_CXXFLAGS="$PKG_CXXFLAGS $MYCXXFLAGS2"
fi

AC_MSG_CHECKING([whether $CXX accepts -ffriend-injection])
AC_LANG([C++])
CXXFLAGS="-ffriend-injection"
AC_COMPILE_IFELSE([[char b[10];]], [ffriend_set=yes], [ffriend_set=no])
AC_MSG_RESULT($ffriend_set)

if test $ffriend_set = yes; then    
   PKG_CXXFLAGS="$PKG_CXXFLAGS -ffriend-injection"
fi

dnl Are we on OSX?
dnl If so, use dmalloc and don't compile our own cblas headers.
AC_MSG_CHECKING([whether we are on Darwin])
UNAME=`uname -a`
echo $UNAME | grep -i Darwin > /dev/null 2>&1
if test "$?" = 0 ; then
 AC_MSG_RESULT(yes)
 cp inst/extras/malloc.c src/malloc.c
 mv src/cblas_dgemm.c inst/extras
 mv src/cblas.h inst/extras
else
 AC_MSG_RESULT(no)
fi

AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_CXXFLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

