#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# != 4 ]; then
  echo "Print the number of common hashes for 2 genomes"
  echo "#1: genome/"
  echo "#2: Genome.id"
  echo "#3: Genome.id"
  echo "#4: CDS|PRT|HMM"
  exit 1
fi
DIR=$1
G1=$2
G2=$3
TYPE=$4


H1=`$THIS/../file2hash $G1`
H2=`$THIS/../file2hash $G2`
F1=$DIR/$H1/$G1/$G1.hash-$TYPE
F2=$DIR/$H2/$G2/$G2.hash-$TYPE
N1=`cat $F1 | wc -l`
N2=`cat $F2 | wc -l`
COMMON=`$THIS/../setIntersect.sh $F1 $F2 1 | wc -l`
echo "$G1 $G2 $N1 $N2 $COMMON"

