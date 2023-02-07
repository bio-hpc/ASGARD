#!/bin/sh

u=`uname`

if [ $u == 'Darwin' ]; then
  (glibtoolize --version) < /dev/null > /dev/null 2>&1 || {
        echo "You must have glibtoolize installed to compile MTK++";
        echo;
        exit;
  }
else
  (libtoolize --version) < /dev/null > /dev/null 2>&1 || {
        echo "You must have libtoolize installed to compile MTK++";
        echo;
        exit;
  }
fi

(automake --version) < /dev/null > /dev/null 2>&1 || {
        echo "You must have automake installed to compile MTK++";
        echo;
        exit;
}

(autoconf --version) < /dev/null > /dev/null 2>&1 || {
        echo "You must have autoconf installed to compile MTK++";
        echo;
        exit;
}

if [ $u == 'Darwin' ]; then
  echo n | glibtoolize --copy --force || exit;
else
  echo n | libtoolize --copy --force || exit;
fi

aclocal || exit;
autoheader || exit;
automake --add-missing --copy;
autoconf || exit;
automake || exit;
#./configure $@

