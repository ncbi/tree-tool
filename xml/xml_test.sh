#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "XML test"
  echo "#1: test viewer (0/1)"
  exit 1
fi
V=$1


cd $THIS

$THIS/xml2schema -qc test/pubmed.xml 
$THIS/xml2schema -qc test/pubmed_30217548.xml 
$THIS/xml2schema -qc test/pubmed.xml 

if [ $V == 1 ]; then
  $THIS/xml_view -qc test/pubmed.xml
fi


TMP=$( mktemp )


$THIS/xml_txt2bin test/pubmed1.xml $TMP.binxml -qc
$THIS/xml_bin2txt $TMP.binxml $TMP.xml -qc
differ $TMP.xml test/pubmed1.xml 


rm $TMP*
