##
## Copyright (C) by Argonne National Laboratory
##     See COPYRIGHT in top-level directory
##

AC_PREREQ([2.69])

m4_include([version.m4])

AC_INIT([MPIFORT], MPICH_VERSION_m4, [discuss@mpich.org], [mpifort], [http://www.mpich.org/])

AC_CONFIG_AUX_DIR(confdb)
AC_CONFIG_MACRO_DIR(confdb)
AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.12.3 subdir-objects no-dist])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_PROG_CC
AM_PROG_CC_C_O

LT_PREREQ([2.2.6])
LT_INIT()

AC_CONFIG_HEADERS([mpifort_config.h])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
