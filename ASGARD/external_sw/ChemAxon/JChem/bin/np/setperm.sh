#!/bin/bash
find ../.. -name "*.sh" -o -name "*.pl" -o -name "*.cgi" -print | xargs chmod a+x
find ../.. -print | grep "\/bin\/" | xargs chmod a+x

