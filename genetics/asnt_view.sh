#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: ASN.1 file"
  exit 1
fi
IN=$1


TMP=$( mktemp )
comment $TMP


$THIS/asnt $IN  -xml $TMP.xml  -qc  1> $TMP
$THIS/../xml/xml_view $TMP.xml -qc


rm $TMP
