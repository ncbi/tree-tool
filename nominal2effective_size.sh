#!/bin/bash
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo -e "${GREEN}Print effective size"
  echo "#1: file with columns: object nominal_value"
  exit 1
fi
IN=$1


S=`cut -f 2 $IN | sort | uniq -c | awk '{s+=$1} END {printf "%f", s}'`
S2=`echo "$S "'*'" $S" | bc -l`
#echo $S2
VAR=`cut -f 2 $IN | sort | uniq -c | awk '{s+=$1*$1/'$S2'} END {printf "%f", s}'`
EFF=`echo 1 / $VAR | bc -l`
echo $EFF
