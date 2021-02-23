#!/bin/bash
if [ $# -ne 1 ]; then
  echo "Exit 1 if there are lineage differences between phen.old/ and phen/"
  echo "#1: Genome.id"
  exit 1
fi
ASM=$1


H=`$THIS/../file2hash $ASM`

diff phen/$H/$ASM phen.old/$H/$ASM | grep '^[><]' | sed 's/^..//1' | sed 's/:.*$//1' | sort | uniq -c | grep -v ' 1 ' 
if [ $? == 0 ]; then
  exit 1
fi



