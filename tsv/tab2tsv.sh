#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Add a header to a tab-delimited file"
  echo "#1: file"
  exit 1
fi
F=$1


if [ ! -s $F ]; then
  exit 0
fi

N=$( head -1 $F | tr '\t' '\n' | wc -l )
i=0
printf "#"
while [ $i -lt $N ]; do
  if [ $i -gt 0 ]; then
    printf '\t'
  fi
  printf "%d" $i
  i=$(( i + 1 ))
done
printf '\n'

cat $F

