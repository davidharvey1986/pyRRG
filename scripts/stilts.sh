#!/bin/sh

#+
#  Name:
#     stilts

#  Purpose:
#     Invokes the STILTS application on unix.

#  Description:
#     This shell script invokes the STILTS application. 
#     It's not very complicated, but performs some argument manipulation
#     prior to invoking java with the right classpath and classname.
#
#     1. if a class path is specified using either the CLASSPATH
#        environment variable or the -classpath flag to this script,
#        it will be added to the application classpath
#
#     2. any initial command-line arguments which look like they are destined
#        for java itself (starting with -D or -X, or prefixed with -J) will
#        be sent to java, and the others will be sent to the application

#  Requisites:
#     - java should be on the path.
#
#     - relative to the directory in which this script is installed,
#       one of the following jar files should exist and contain the
#       STILTS classes:
#          ./stilts.jar
#          ../../lib/ttools/stilts-app.jar
#          ./stilts-app.jar
#          ./topcat-full.jar
#          ./topcat-lite.jar
#       (on a Mac it looks in the resource bundle expected from a dmg
#       installation as well)

#  Authors:
#     MBT: Mark Taylor (Starlink)
#-

#  Find where this script is located.
scriptname="$0"
while [ -L "$scriptname" ]
do
   scriptname=`readlink "$scriptname" 2>/dev/null`
done
test -n "$scriptname" || scriptname=$0
bindir="`dirname $scriptname`"

#  Set locations of acceptable jar files (relative to this script).
lib_path=`pip show pyRRG | grep Location | cut -f2 -d' '`
stilts_jars="\
 $lib_path/stilts/stilts.jar\
 $bindir/stilts.jar\
 $bindir/../lib/ttools/stilts-app.jar\
 $bindir/stilts-app.jar\
 $bindir/topcat-full.jar\
 $bindir/topcat-lite.jar\
"

# If we're on a Mac, look for the jar in the resource bundle,
# since we may have been installed that way.  This is a bit of a mess,
# because of historical changes in the locations of these scripts.
if test -x /usr/bin/sw_vers && /usr/bin/sw_vers | grep -iq 'Mac *OS'; then
   bundle_jar="Java/topcat-full.jar"
   stilts_jars="$stilts_jars \
                $bindir/../../$bundle_jar"

fi

#  Locate the application jar file.
for j in $stilts_jars; do
   if test -z "$appjar" -a -f "$j"; then
      appjar="$j"
   fi
done
if test ! -f "$appjar"
then
   echo 1>&2 "Can't find stilts classes relative to ${scriptname} - looked for:"
   echo " " "$stilts_jars" | sed 1>&2 's/  */\n   /g'
   exit 1
fi

#  Pull out any arguments which look to be destined for the java binary.
javaArgs=""
while test "$1"
do
   if echo $1 | grep -- '^-[XD]' >/dev/null; then
      javaArgs="$javaArgs $1"
      shift
   elif echo $1 | grep -- '^-J' >/dev/null; then
      javaArgs="$javaArgs `echo $1 | sed s/^-J//`"
      shift
   elif [ "$1" = "-classpath" -a -n "$2" ]; then
      shift
      export CLASSPATH="$1"
      shift
   else
      break
   fi
done

#  Check for Cygwin and transform paths.
case "`uname`" in
  CYGWIN*)
    if test -n "$CLASSPATH"; then
       CLASSPATH=`cygpath --path --windows "${appjar}:$CLASSPATH"`
    else
       CLASSPATH=`cygpath --windows "${appjar}"`
    fi
  ;;
  *)
    CLASSPATH="${appjar}:${CLASSPATH}"
  ;;
esac

# Execute the command.
java \
   $javaArgs \
   -classpath $CLASSPATH \
   uk.ac.starlink.ttools.Stilts \
   "$@"
