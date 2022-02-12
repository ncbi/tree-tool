#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "tsv_join: #1 + #2 -> #1"
  echo "create #1.bak"
  echo "#1: tsv-file name 1"
  echo "#2: tsv-file name 2"
  echo "#3: tsv_joinoption: -left|-remove|''"
  exit 1
fi
F1=$1
F2=$2
OPT="$3"


mv $F1 $F1.bak
$THIS/tsv_join $F1.bak $F2 -qc $OPT > $F1



