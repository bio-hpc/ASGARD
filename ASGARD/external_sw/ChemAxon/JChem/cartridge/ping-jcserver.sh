#/bin/sh

#####################################################################
# Checks if a JChem (Cartridge) Server is available on the specified
# host at the specified port.
#
# Exit codes:
#    0: a JChem Server is available and is of the same version as the
#       jchem.jar pointed to by the "jchem_jar" variable.
#    1: a JChem Server is available and is of a version different 
#       from the jchem.jar pointed to by the "jchem_jar" variable.
#    2: no service is available on the specified host at the specified
#       port.
#    3: an unknown service is listening on the specified host at the
#       specified port.
#    4: an undefined error occurred while trying to ping the service.
#
# Usage:
#    ping-jcserver.sh <hostname>:<port>
#
# Example:
#    ping-jcserver.sh localhost:1099
#
#####################################################################

jchem_jar=../lib/jchem.jar

"$JAVA_HOME/bin/java" \
    -classpath $jchem_jar \
    -Djava.util.logging.config.class=chemaxon.jchem.cartridge.util.LoggingConfigurator \
    -Dchemaxon.jchem.cartridge.config.file=conf/jcart.properties \
    chemaxon.jchem.cartridge.install.JChemServerPinger $@

