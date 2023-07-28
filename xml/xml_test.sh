#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "XML test"
  echo "#1: go"
  exit 1
fi


cd $THIS

$THIS/xml2schema -qc pubmed.xml 
$THIS/xml2schema -qc pubmed_30217548.xml 
$THIS/xml2schema -qc fig9.xml 

$THIS/xml_view -qc fig9.xml


TMP=`mktemp`

$THIS/xml_txt2bin fig9.xml $TMP.binxml -qc
$THIS/xml_bin2txt $TMP.binxml $TMP.xml -qc
diff fig9.xml $TMP.xml

rm $TMP*
