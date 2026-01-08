#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# != 3 ]; then
  echo "Print common hashes with frequences"
  echo "#1: genome/"
  echo "#2: List of Genome.id's"
  echo "#3: CDS|PRT|HMM"
  exit 1
fi
DIR=$1
LIST=$2
TYPE=$3


$THIS/../../trav $LIST -step 1 "cat $DIR/%h/%f/%f.hash-$TYPE" | sort | uniq -c | sort -k1nr,1 -k2n,2


