#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# != 5 ]; then
  echo "Print the number of common hashes for 2 genomes"
  echo "#1: genome/"
  echo "#2: Genome.id"
  echo "#3: Genome.id"
  echo "#4: CDS|PRT|HMM"
  echo "#5: out file with common hashes or ''"
  exit 1
fi
DIR=$1
G1=$2
G2=$3
TYPE=$4
OUT="$5"


H1=$( $THIS/../../file2hash $G1 )
H2=$( $THIS/../../file2hash $G2 )
F1=$DIR/$H1/$G1/$G1.hash-$TYPE
F2=$DIR/$H2/$G2/$G2.hash-$TYPE
N1=$( < $F1  wc -l )
N2=$( < $F2  wc -l )
COMMON=$( $THIS/../../setIntersect.sh $F1 $F2 1 | wc -l )
echo -e "$G1\t$G2\t$N1\t$N2\t$COMMON"

if [ "$OUT" ]; then
  $THIS/../../setIntersect.sh $F1 $F2 1 > $OUT
fi

