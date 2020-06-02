#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print: #1 trunc_n trunc_c"
  echo "#1: Protein accession"
  exit 1
fi
ACC=$1


TMP=`mktemp`
#echo $TMP  


wget 'https://www.ncbi.nlm.nih.gov/Structure/cdannots/cdannots.fcgi?queries='$ACC  -O $TMP.html  -o $TMP.wget
$THIS/../cddJson $TMP.html > $TMP.trunc
L=(`cat $TMP.trunc`)

echo "$ACC ${L[*]}"


rm $TMP*  
