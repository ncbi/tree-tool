#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print a .dm-file"
  echo "#1: File with numeric columns"
  echo "#2: Number precision"
  echo "#3: Named (0/1)"
  exit 1
fi


N=`cat $1 | wc -l`
if [ $N == 0 ]; then
  echo "No data"
  exit 1
fi

C=`head -1 $1 | wc -w`
if [ $3 == 1 ]; then
  NAME="name"
  C=$(( $C - 1 ))
else
  NAME="noname"
fi

echo "ObjNum $N $NAME nomult" 
echo "Attributes"
i=1
while [ $i -le $C ]; do
  echo "  V$i real $2"
  i=$(( $i + 1 ))
done
echo "Data"
cat $1

