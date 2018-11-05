dnl $1 is src/mpid/ch4/netmod/ucx/ucx
AC_DEFUN([PAC_SUBDIR_UCX],[
    PAC_SET_HEADER_LIB_PATH(ucx)
    ucx_embedded=""
    if test "${with_ucx}" = "embedded" ; then
        ucx_embedded="yes"
    elif test -z ${with_ucx} ; then
        if test -f "$1/configure" ; then
            ucx_embedded="yes"
        else
            ucx_embedded="no"
        fi
    else
        ucx_embedded="no"
    fi

    if test "${ucx_embedded}" = "yes" ; then
        PAC_SUBDIR_EMBED([ucp],[$1],[$1/src/ucp],[$1/src])
        PAC_CONFIG_ARG_UCX
    else
        PAC_PUSH_FLAG(LIBS)
        PAC_CHECK_HEADER_LIB_FATAL(ucx, ucp/api/ucp.h, ucp, ucp_config_read)
        PAC_POP_FLAG(LIBS)
        PAC_SUBDIR_SYSTEM([ucp])
        PAC_SUBDIR_SYSTEM([ucs])
    fi
])

AC_DEFUN([PAC_CONFIG_ARG_UCX],[
    ucx_config_arg="--disable-static --enable-embedded"
    export ucx_config_arg
])
