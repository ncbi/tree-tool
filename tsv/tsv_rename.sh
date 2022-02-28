#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Rename a column in place, create #1.bak"
  echo "#1: tsv-file name"
  echo "#2: column number (>=1)"
  echo "#3: new column name"
  exit 1
fi
F=$1
N=$2
NAME="$3"


mv $F $F.bak
$THIS/tsv_rename -qc $F.bak $N "$NAME" > $F



