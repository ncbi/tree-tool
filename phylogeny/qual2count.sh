#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo -e "Print pairs: taxon\tgains+losses"
  echo "cut -f 2 | count: sum - count = misconguence"
  echo "#1: qual-file made by makeFeatureTree"
  exit 1
fi
QUAL=$1


cat $QUAL | sed 's| / .*$||1' | sed 's/ [-+]/\t/g' | tr ' ' '_' | awk -F '\t' '{printf "%s\t%d\n", $1, $2 + $3};' 


