#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Print a .dm-file with a two-way attribute 'dist'"
  echo "#1: File with distances in a Phylip format"
  echo "#2: max. number of decimals in distances"
  exit 1
fi
F=$1
DEC=$2


L=`cat $F | wc -l`
if [ $L == 0 ]; then
  error "No data"
fi

N=`head -1 $F`
M=$(( $L - 1 ))
if [ $N -ne $M ]; then
  error "Number of objects does not match the number of lines"
fi

echo "ObjNum $N name nomult" 
echo "Attributes"
echo "  dist positive2 $DEC"
echo "DATA"
tail -n +2 $F | cut -f 1 -d ' '
echo "dist FULL"
tail -n +2 $F | cut -f 1 -d ' ' --complement 
