# Check for QT compiler flags, linker flags, and binary packages
# This was adapted from autoqt for QT4+
AC_DEFUN([CHECK_QT],
[
AC_REQUIRE([AC_PROG_CXX])
AC_REQUIRE([AC_PATH_X])

AC_MSG_CHECKING([QTDIR])
AC_ARG_WITH([qtdir], [  --with-qtdir=DIR        Qt installation directory [default=$QTDIR]], QTDIR=$withval)

# Check that QTDIR is defined or that --with-qtdir given
#if test x"$QTDIR" = x ; then
#
#    QT_SEARCH="/usr/local/Trolltech/Qt-4.3.3"
#    for i in $QT_SEARCH; do
#        if test -f $i/include/Qt/qglobal.h -a x$QTDIR = x; then QTDIR=$i; fi
#    done
#fi

if test x"$QTDIR" = x ; then
    AC_MSG_WARN([QTDIR is not defined])
    test_qt="no"
else
    AC_MSG_RESULT([$QTDIR])
    test_qt="yes"

    # Change backslashes in QTDIR to forward slashes to prevent escaping
    # problems later on in the build process, mainly for Cygwin build
    # environment using MSVC as the compiler
    QTDIR=`echo $QTDIR | perl -p -e 's/\\\\/\\//g'`

    # Figure out which version of Qt we are using
    AC_MSG_CHECKING([Qt version])

    QT_VER=`grep 'define.*QT_VERSION_STR\W' $QTDIR/include/Qt/qglobal.h | perl -p -e 's/\D//g'`
    case "${QT_VER}" in
      4*)
         QT_MAJOR="4"
      ;;
      *)
         AC_MSG_ERROR([*** Don't know how to handle this Qt major version])
      ;;
    esac
    AC_MSG_RESULT([$QT_VER ($QT_MAJOR)])

    # Check that moc is in path
    AC_CHECK_PROG(MOC, moc, moc)
    if test x$MOC = x ; then
        AC_MSG_WARN([moc not found])
    fi

    # uic is the Qt user interface compiler
    AC_CHECK_PROG(UIC, uic, uic)
    if test x$UIC = x ; then
        AC_MSG_WARN([uic not found])
    fi

    # On unix, figure out if we're doing a static or dynamic link
    case "${host}" in
      *apple*)
        QT_IS_STATIC=`ls $QTDIR/lib/*.la 2> /dev/null`
        if test "x$QT_IS_STATIC" = x; then
          QT_IS_STATIC="yes"
        else
          QT_IS_STATIC="no"
        fi
        if test x$QT_IS_STATIC = xno ; then
          QT_IS_DYNAMIC=`ls $QTDIR/lib/QtCore.framework/QtCore 2> /dev/null`
          if test "x$QT_IS_DYNAMIC" = x;  then
            AC_MSG_WARN([Couldn't find any Qt libraries])
          else
            QT_IS_DYNAMIC="yes"
          fi
        fi
        ;;
      *-*-linux*)
        QT_IS_STATIC=`ls $QTDIR/lib/libQtCore*.a 2> /dev/null`
        if test "x$QT_IS_STATIC" = x; then
          QT_IS_STATIC="no"
        else
          QT_IS_STATIC="yes"
        fi
        if test x$QT_IS_STATIC = xno ; then
          QT_IS_DYNAMIC=`ls $QTDIR/lib/*.so 2> /dev/null`
          if test "x$QT_IS_DYNAMIC" = x;  then
            AC_MSG_WARN([Couldn't find any Qt libraries])
          else
            QT_IS_DYNAMIC="yes"
          fi
        fi
        ;;
    esac

    # Calculate Qt include path
    QT_CXXFLAGS="-I$QTDIR/include"
    QT_DEFINES="-DUSE_QT -DQT_CORE_LIB -DQT_XML_LIB -DQT_SHARED "
    QT_LIBS=""
    QT_INCPATH=""

    case "${host}" in
      *apple*)
        if test $QT_IS_DYNAMIC = yes ; then
            QT_INCPATH="-I$QTDIR/mkspecs/macx-g++ -I. -I$QTDIR/lib/QtCore.framework/Versions/4/Headers -I$QTDIR/include/QtCore -I$QTDIR/include/QtCore -I$QTDIR/lib/QtXml.framework/Versions/4/Headers -I$QTDIR/include/QtXml -I$QTDIR/include/QtXml -I$QTDIR/include -I. -Imoc -I. -F$QTDIR/lib"
            QT_LIBS="-F$QTDIR/lib -L$QTDIR/lib"
            QT_LFLAGS="-headerpad_max_install_names"
        fi
        ;;
    *-*-linux*)
        if test $QT_IS_DYNAMIC = yes ; then
            QT_INCPATH="-I$QTDIR/mkspecs/linux-g++ -I. -I$QTDIR/include/QtCore -I$QTDIR/include/QtXml -I$QTDIR/include -I. -Imoc -I."
            QT_LIBS="-L$QTDIR/lib -lQtCore -lQtXml"
            QT_LFLAGS="-headerpad_max_install_names"
        fi
        ;;
    esac

    AC_MSG_CHECKING([QT_CXXFLAGS])
    AC_MSG_RESULT([$QT_CXXFLAGS])
    AC_MSG_CHECKING([QT_DEFINES])
    AC_MSG_RESULT([$QT_DEFINES])
    AC_MSG_CHECKING([QT_LIBS])
    AC_MSG_RESULT([$QT_LIBS])
    AC_MSG_CHECKING([QT_INCPATH])
    AC_MSG_RESULT([$QT_INCPATH])

    AC_SUBST(QT_CXXFLAGS)
    AC_SUBST(QT_DEFINES)
    AC_SUBST(QT_LIBS)
    AC_SUBST(QT_INCPATH)
fi
])
