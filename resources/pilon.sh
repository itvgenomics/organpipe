#! /bin/sh
set -e

# export JAVA_HOME=/usr/lib/jvm/default-java
export JAVA_CMD=java

# Include the wrappers utility script
. /usr/lib/java-wrappers/java-wrappers.sh

# For memory setting see https://github.com/rrwick/Unicycler/issues/63
run_java -Xmx50G -Xss2560k -jar /usr/share/java/pilon.jar "$@"
