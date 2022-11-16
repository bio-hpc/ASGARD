#!/bin/bash
CLASSPATH=`echo "$CLASSPATH"|sed 's/jchem\.jar/nojchem.jar/g'`
export CLASSPATH
exec "$@"

