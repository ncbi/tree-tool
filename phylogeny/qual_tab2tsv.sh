#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print qual-file in .tsv format"
  echo "#1: qual tab-file"
  exit 1
fi
TAB=$1


echo -e "#Feature\tgains\tlosses\tcriterion\tobjects\toptional" 
awk -F '\t' '{OFS="\t"; print $1, $2, $3, $2 + $3 - 1, $4, $5};' $TAB

