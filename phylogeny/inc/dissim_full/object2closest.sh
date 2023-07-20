#!/bin/bash --noprofile
THIS=`dirname $0`
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print 100 approximately closest objects to #1"
  echo "#1: object id"
  echo "#2: object directory or '' "
  echo "Requires: Time: O(log^4(n))"
  exit 1
fi
OBJ=$1
DIR="$2"


if [ $DIR ]; then
  error "Using $DIR is not implemented"
fi


TMP=`mktemp`
#echo $TMP
#set -x


# $TMP.other
awk -F '-' '{if ($1 == "'$OBJ'") print $2, $3};' $THIS/../dissim_full >  $TMP.select
awk -F '-' '{if ($2 == "'$OBJ'") print $1, $3};' $THIS/../dissim_full >> $TMP.select
sort $TMP.select > $TMP.other

CPP_DIR/phylogeny/tree2obj.sh $THIS/tree > $TMP.list
join  -1 1  -2 1  $TMP.list $TMP.other | sort -k 2 -g > $TMP.sort 
head -100 $TMP.sort | cut  -f 1  -d ' ' 


rm $TMP*
