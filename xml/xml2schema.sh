#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print derived XML schema"
  echo "#1: gzipped XML file"
  exit 1
fi
BIOS=$1


TMP=`mktemp`

gunzip $BIOS -c > $TMP
$THIS/xml2schema $TMP -qc 

rm $TMP*
