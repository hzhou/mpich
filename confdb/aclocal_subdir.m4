dnl This is parallel to aclocal_subcfg.m4.
dnl `aclocal_subcfg.m4` executes subcfg inside top configure's space.
dnl `aclocal_subdir.m4`(this file) runs separate `configure` scripts.
dnl Both files have distinct logic despite similar sematics, therefore separate m4.
dnl FIXME: PAC_CONFIG_SUBDIR_ARGS and PAC_CONFIG_SUBDIR need be moved here.
dnl        Also: PAC_SUBDIR_CACHE_CLEANUP

dnl In MPICH, there are: 
dnl     src/mpl
dnl     src/openpa
dnl     src/izem
dnl     src/mpi/romio
dnl     src/hwloc
dnl     [devicedir] # none left, correct?
dnl     src/pm/[pm_name]
dnl     src/nameserv/[namepublisher]
dnl     src/util/logging/[logging_name]
dnl     src/*/libfabric
dnl     src/*/ucx
dnl     what about bindingsubsystems

dnl call this right after AC_INIT (not must, but a good convention).
AC_DEFUN([PAC_SUBDIR_INIT],[
    dnl automake will define all AC_SUBST's
    AC_SUBST([external_subdirs])
    AC_SUBST([external_dist_subdirs])
    AC_SUBST([external_libs])
])

dnl call this at the end of configure.ac
dnl This will inflate the configure for hydra and romio, unfortunately
AC_DEFUN([PAC_CONFIG_ALL_SUBDIRS],[
    all="m4_default([$1], [$subsystems])"
    echo ********************************************
    echo Configure all subsystems: $all
    echo ********************************************
    if test -n "$all" ; then
        for subsys in $all ; do 
            case $subsys in
                mpl|*/mpl)
                    config_arg="$mpl_config_arg"
                    ;;
                */openpa)
                    config_arg="$opa_config_arg"
                    ;;
                */izem)
                    config_arg="$zm_config_arg"
                    ;;
                */hwloc)
                    config_arg="$hwloc_config_arg"
                    ;;
                */libfabric)
                    config_arg="$ofi_config_arg"
                    ;;
                */ucx)
                    config_arg="$ucx_config_arg"
                    ;;
                src/mpi/romio)
                    config_arg="--with-mpl-prefix=../../mpl"
                    ;;
                src/pm/hydra)
                    config_arg="--with-mpl-prefix=../../mpl --with-hwloc-prefix=../../hwloc"
                    ;;
                *)
                    config_arg=""
                    ;;
            esac

            dnl reset all flags, if specific package need inherit config time flags, do so explicitly
            PAC_PUSH_ALL_FLAGS
            PAC_RESET_ALL_FLAGS
            PAC_CONFIG_SUBDIR_ARGS([$subsys],[$config_arg],[],[AC_MSG_ERROR([$subsys configure failed])])
            PAC_POP_ALL_FLAGS
        done 
    fi

    dnl FIXME: Not sure why we need following. Comment out to find out the consequence...
    dnl if test "$DEBUG_SUBDIR_CACHE" = yes -a "$enable_echo" != yes ; then 
    dnl     set +x
    dnl fi
    dnl PAC_SUBDIR_CACHE_CLEANUP

    dnl # Make subsystems available to makefiles.
    dnl # FIXME does the makefile actually need this?
    dnl subsystems="$devsubsystems $subsystems $bindingsubsystems"
])

AC_DEFUN([PAC_CONFIG_ALL_PARALLEL], [
    all="m4_default([$1], [$subsystems])"
    romio_config_arg="--with-mpl-prefix=../../mpl"
    hydra_config_arg="--with-mpl-prefix=../../mpl --with-hwloc-prefix=../../hwloc"
    export romio_config_arg hydra_config_arg
    j=m4_default([$2], 2)
    perl maint/subdir_config.pl -j$j $all
])

dnl ****************************************
dnl Apparently we need allow the possiblity of library under custom name.
dnl PAC_SUBDIR_LIBNAME([var], [default])
AC_DEFUN([PAC_SUBDIR_LIBNAME],[
    AC_ARG_VAR([$1],[can be used to override the name of the $2 library (default: "$2")])
    $1=${$1:-"$2"}
    export $1
    AC_SUBST($1)
])

dnl There are three distinct possibilities: embedded, system, user supplied prefix.

dnl **** embed ************************************
dnl PAC_SUBDIR_EMBED([libname],[src_dir],[la_dir],[inc_dir])
AC_DEFUN([PAC_SUBDIR_EMBED],[
    if test -e "$2" ; then
        subsystems="$subsystems $2"
        external_subdirs="$external_subdirs $2"
        external_dist_subdirs="$external_dist_subdirs $2"
        external_libs="$external_libs $3/lib$1.la"
        PAC_APPEND_FLAG([-I./$4],[CPPFLAGS])
        if test $srcdir != . ; then
            PAC_APPEND_FLAG([-I${srcdir}/$4],[CPPFLAGS])
        fi
    else
        dnl FIXME: if may fail, then should fail here.
        AC_MSG_WARN([Attempted to use the embedded $1 source tree in "$2", but it is missing.  Configuration or compilation may fail later.])
    fi
])

dnl Reuse the convenience library, skip subsystems
dnl PAC_SUBDIR_EMBED([libname],[src_dir],[la_dir],[inc_dir])
dnl For sub library such as ROMIO, it may not want to link with certain library, 
dnl such as mpl, or it may end up linking it twice in the upper layer. 
dnl `m4_define([EMBED_NO_LIB], 1) to skip it. Remember to undefine it right away.
dnl
AC_DEFUN([PAC_SUBDIR_EMBED_REUSE],[
    if test -e "$2" ; then
        m4_ifndef([EMBED_NO_LIB],[
            external_libs="$external_libs $3/lib$1.la"
        ],
        []
        )
        PAC_APPEND_FLAG([-I./$4],[CPPFLAGS])
        if test $srcdir != . ; then
            PAC_APPEND_FLAG([-I${srcdir}/$4],[CPPFLAGS])
        fi
    else
        dnl FIXME: if it may fail later, then it should fail here.
        AC_MSG_WARN([Attempted to reuse the embedded $1 source tree in "$2", but it is missing.  Configuration or compilation may fail later.])
    fi
])

dnl **** system ************************************
dnl PAC_SUBDIR_SYSTEM(libname)
AC_DEFUN([PAC_SUBDIR_SYSTEM],[
    PAC_PREPEND_FLAG([-l$1],[LIBS])
    m4_ifset([INSIDE_MPICH],[
            PAC_PREPEND_FLAG([-l$1],[WRAPPER_LIBS])
        ],[])
])

dnl **** prefix ************************************
dnl PAC_ARG_WITH_PREFIX(libname)
AC_DEFUN([PAC_ARG_WITH_PREFIX],[
	AC_ARG_WITH([$1-prefix],
            [AS_HELP_STRING([[--with-$1-prefix[=DIR]]], [use the $1
                            library installed in DIR, rather than the
                            one included in the distribution.  Pass
                            "embedded" to force usage of the included
                            $1 source.])],
            [],dnl action-if-given
            [with_$1_prefix="embedded"])
])

dnl PAC_SUBDIR_PREFIX([libname],[prefix_dir])
AC_DEFUN([PAC_SUBDIR_PREFIX],[
    PAC_APPEND_FLAG([-I$2/include],[CPPFLAGS])
    if test -d "$2/lib64"; then
        pac_libdir="$2/lib64"
    else
        pac_libdir="$2/lib"
    fi
    PAC_PREPEND_FLAG([-L${pac_libdir}],[LDFLAGS])
    m4_ifset([INSIDE_MPICH],[
        PAC_PREPEND_FLAG([-L${pac_libdir}],[WRAPPER_LDFLAGS])
        ],[])

    if "$3"x == x ; then
        PAC_SUBDIR_SYSTEM([$1])
    fi
])

dnl .internal.
dnl PAC_SUBDIR_PREFIX_CHECK([name],[dir],[file])
AC_DEFUN([PAC_SUBDIR_PREFIX_CHECK],[
    AS_IF([test -s "$2/$3"],
        [:],[AC_MSG_ERROR([the $1 installation in "$2" appears broken])])
])

