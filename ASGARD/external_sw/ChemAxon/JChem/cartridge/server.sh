#!/bin/sh

java_debug=false
trace_lichandler=false
other=

if [ -z "$JAVA_HOME" ];
then
    echo "JAVA_HOME is not set";
    exit 1;
elif [ ! -x "$JAVA_HOME"/bin/java ];
then
    echo "Invalid JAVA_HOME environment variable" 
    echo "No executable $JAVA_HOME/bin/java found"
    exit 1;
else
    jcsrv_jvm="$JAVA_HOME"/bin/java;
fi

[ "$java_debug" = "true" ] && debug=-agentlib:jdwp=transport=dt_socket,address=30303,server=y,suspend=n
[ "$trace_lichandler" = "true" ] && lichandtrace='-Dchemaxon.license.verify.verbose=cfgfile'
other_options="$debug $lichandtrace -XX:-OmitStackTraceInFastThrow"

if [ "$OSTYPE" = "cygwin" ];
then
    PATHSEP=';'
    npa() {
        cygpath -wp $1
    }
else
    PATHSEP=':'
    npa() {
        echo $1;
    }
fi;

cp=$(dirname $0)/../lib/jchem.jar


configdir=$(dirname $0)/conf

set -x
"$jcsrv_jvm" \
   -Djava.util.logging.config.class=chemaxon.jchem.cartridge.util.LoggingConfigurator \
   -Dchemaxon.jchem.cartridge.config.file=conf/jcart.properties \
   -Djava.awt.headless=true \
   -classpath $cp \
   $other_options \
   chemaxon.jchem.cartridge.server.Bootstrapper $@
set +x
