dnl $1 is src/mpid/ch4/netmod/ofi/libfabric
AC_DEFUN([PAC_SUBDIR_OFI],[
    PAC_SET_HEADER_LIB_PATH(libfabric)
    ofi_embedded=""
    AC_SUBST([ofi_embedded])
    AC_SUBST([with_libfabric])
    if test "${with_libfabric}" = "embedded" ; then
        ofi_embedded="yes"
    elif test -z "${with_libfabric}" && test -z "${with_libfabric_lib}" && test -z "${with_libfabric_include}" ; then
        if test -f $1/configure ; then
            ofi_embedded="yes"
        else
            ofi_embedded="no"
        fi
    else
        ofi_embedded="no"
    fi

    if test "{ofi_embedded}" = "yes" ; then
        AC_MSG_NOTICE([CH4 OFI Netmod:  Using an embedded libfabric])
        PAC_SUBDIR_EMBED([libfabric],[$1],[$1/src],[$1/include])
        if test x"$ofi_direct_provider" != x ; then
            PAC_APPEND_FLAG([-I${builddir}/$1/prov/${ofi_direct_provider}/include],[CPPFLAGS])
            PAC_APPEND_FLAG([-I${srcdir}/$1/prov/${ofi_direct_provider}/include],[CPPFLAGS])
            PAC_APPEND_FLAG([-DFABRIC_DIRECT],[CPPFLAGS])
        fi
        PAC_CONFIG_ARG_OFI
    else
        AC_MSG_NOTICE([CH4 OFI Netmod:  Using an external libfabric])
        PAC_CHECK_SYSTEM_OFI
        PAC_SUBDIR_SYSTEM([fabric])
    fi
])

AC_DEFUN([PAC_CHECK_SYSTEM_OFI],[
    PAC_PUSH_FLAG(LIBS)
    PAC_CHECK_HEADER_LIB([rdma/fabric.h], [fabric], [fi_getinfo], [have_libfabric=yes], [have_libfabric=no])
    PAC_POP_FLAG(LIBS)
    if test "$have_libfabric" = "no" ; then
        AC_MSG_ERROR([Provided libfabric installation (--with-libfabric=${with_libfabric}) could not be configured.])
    fi
])
AC_DEFUN([PAC_CONFIG_ARG_OFI],[
    dnl AC_MSG_NOTICE([RUNNING CONFIGURE FOR ch4:ofi])
    ofi_config_arg="--enable-embedded"
    # set --enable-xxx based on $netmod_args (ref: src/mpid/ch4/subconfigure.m4)
    PAC_OFI_RUNTIME_CAPABILITIES
    if test x"${netmod_args}" != x ; then
        ofi_config_arg="$ofi_config_arg --enable-psm=${enable_psm}"
        ofi_config_arg="$ofi_config_arg --enable-psm2=${enable_psm2}"
        ofi_config_arg="$ofi_config_arg --enable-sockets=${enable_sockets}"
        ofi_config_arg="$ofi_config_arg --enable-verbs=${enable_verbs}"
        ofi_config_arg="$ofi_config_arg --enable-usnic=${enable_usnic}"
        ofi_config_arg="$ofi_config_arg --enable-gni=${enable_gni}"
        ofi_config_arg="$ofi_config_arg --enable-bgq=${enable_bgq}"
        ofi_config_arg="$ofi_config_arg --enable-udp=${enable_udp}"
        ofi_config_arg="$ofi_config_arg --enable-rxm=${enable_rxm}"
        ofi_config_arg="$ofi_config_arg --enable-rxd=${enable_rxd}"
        ofi_config_arg="$ofi_config_arg --enable-tcp=${enable_tcp}"
        ofi_config_arg="$ofi_config_arg --enable-shm=${enable_shm}"
        ofi_config_arg="$ofi_config_arg --enable-mlx=${enable_mlx}"
    fi
    if test "x${ofi_direct_provider}" != "x" ; then
        AC_MSG_NOTICE([Enabling direct embedded provider: ${ofi_direct_provider}])
        ofi_config_arg="$ofi_config_arg --enable-direct=${ofi_direct_provider}"
    fi
    export ofi_config_arg
])

dnl .internal.
dnl $netmod_args - contains the OFI provider
AC_DEFUN([PAC_OFI_RUNTIME_CAPABILITIES],[
    runtime_capabilities="no"
    no_providers="no"
    if test "${netmod_args#* }" != "$netmod_args" ; then
        runtime_capabilities="yes"
        AC_MSG_NOTICE([Using runtime capability set due to multiple selected providers])
    elif test "x$netmod_args" = "x" || test "$netmod_args" = "runtime"; then
        runtime_capabilities="yes"
        no_providers="yes"
        AC_MSG_NOTICE([Using runtime capability set due to no selected provider or explicity runtime selection])
    fi

    if test "$no_providers" = "no" ; then
        enable_psm="no"
        enable_psm2="no"
        enable_sockets="no"
        enable_verbs="no"
        enable_usnic="no"
        enable_gni="no"
        enable_bgq="no"
        enable_udp="no"
        enable_rxm="no"
        enable_rxd="no"
        enable_tcp="no"
        enable_shm="no"
        enable_mlx="no"
    else
        enable_psm="yes"
        enable_psm2="yes"
        enable_sockets="yes"
        enable_verbs="yes"
        enable_usnic="yes"
        enable_gni="yes"
        enable_bgq="yes"
        enable_udp="yes"
        enable_rxm="yes"
        enable_rxd="yes"
        enable_tcp="yes"
        enable_shm="yes"
        enable_mlx="yes"
    fi

    for provider in $netmod_args ; do
        case "$provider" in
            dnl For these providers, we know which capabilities we want to select by default.
            dnl We want to enable compiling with them, but we don't need to enable runtime checks.
            "psm")
                enable_psm="yes"
                ;;
            "psm2" | "opa")
                enable_psm2="yes"
                ;;
            "sockets")
                enable_sockets="yes"
                ;;
            "gni")
                enable_gni="yes"
                ;;
            "bgq")
                enable_bgq="yes"
                ;;

            dnl For these providers, we don't know exactly which capabilities we
            dnl want to select by default so we turn on runtime checks. At some point
            dnl in the future, we may create a specific capability set for them.
            "verbs")
                enable_verbs="yes"
                runtime_capabilities="yes"
                ;;
            "usnic")
                enable_usnic="yes"
                runtime_capabilities="yes"
                ;;
            "udp")
                enable_udp="yes"
                runtime_capabilities="yes"
                ;;
            "rxm")
                enable_rxm="yes"
                runtime_capabilities="yes"
                ;;
            "rxd")
                enable_rxd="yes"
                runtime_capabilities="yes"
                ;;
            "tcp")
                enable_tcp="yes"
                runtime_capabilities="yes"
                ;;
            "shm")
                enable_shm="yes"
                runtime_capabilities="yes"
                ;;
            "mlx")
                enable_mlx="yes"
                runtime_capabilities="yes"
                ;;
            *)
                AC_MSG_WARN("Invalid provider $provider")
        esac
    done

    if test "$runtime_capabilities" = "yes" ; then
        AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
    else
        case "$netmod_args" in
            "psm")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_PSM], [1], [Define to use PSM capability set])
                enable_psm="yes"
                ;;
            "psm2" | "opa")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_PSM2], [1], [Define to use PSM2 capability set])
                enable_psm2="yes"
                ;;
            "sockets")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_SOCKETS], [1], [Define to use sockets capability set])
                enable_sockets="yes"
                ;;
            "gni")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_GNI], [1], [Define to use gni capability set])
                enable_gni="yes"
                ;;
            "bgq")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_BGQ], [1], [Define to use bgq capability set])
                enable_bgq="yes"
                ;;
            "verbs")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_verbs="yes"
                ;;
            "usnic")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_usnic="yes"
                ;;
            "udp")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_udp="yes"
                ;;
            "rxm")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_rxm="yes"
                ;;
            "rxd")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_rxd="yes"
                ;;
            "tcp")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_tcp="yes"
                ;;
            "shm")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_shm="yes"
                ;;
            "mlx")
                AC_DEFINE([MPIDI_CH4_OFI_USE_SET_RUNTIME], [1], [Define to use runtime capability set])
                enable_mlx="yes"
                ;;
            *)
                AC_MSG_WARN("Invalid provider $netmod_args")
        esac
    fi

])

