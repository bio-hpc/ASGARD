#!/bin/bash

JAVA_CMD="$JAVA_HOME/bin/java"
if [ ! -x "$JAVA_CMD" ];
then
    echo "JAVA_HOME is not appropriately set";
    exit 1;
fi

if [ ! -f ../lib/jchem.jar ];
then
    echo "You current working directory is not appropriate."
    exit 1;
fi

set -x
$JAVA_CMD \
    -classpath ../lib/jchem.jar \
    -Dchemaxon.jchem.cartridge.config.file=conf/jcart.properties \
    -Djava.util.logging.config.class=chemaxon.jchem.cartridge.util.LoggingConfigurator \
    -Dchemaxon.jchem.cartridge.logging.config.file=conf/upgrade-logging.properties \
    chemaxon.jchem.cartridge.install.UpgradeCmdLine $@
set +x
