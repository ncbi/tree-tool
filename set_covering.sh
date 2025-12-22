#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Set-covering problem (approximation)"
  echo "Print a minimal covering set"
  echo "#1: list of sets: each row is a space-delimited list (set)"
  exit 1
fi
IN=$1


TMP=$( mktemp )
#comment $TMP


mkdir $TMP.sets
i=1
while read -r SET
do
  echo $SET | tr ' ' '\n' | grep -vx '' | sort -u > $TMP.sets/$i || true
  i=$(( i + 1 ))
done < $IN

while true
do
  $THIS/trav -noprogress $TMP.sets "cat %d/%f" > $TMP.all
  if [ ! -s $TMP.all ]; then
    break
  fi
  sort $TMP.all | uniq -c | sort -k1nr | head -1 > $TMP.1
  L=( $( cat $TMP.1 ) )
  N=${L[0]}
  W=${L[1]}
  if [ $N == 1 ]; then
    break;
  fi
  printf "\r%d " $N > /dev/stderr
  echo $W
  ls $TMP.sets > $TMP.list
  $THIS/trav -noprogress $TMP.list "grep -xq $W $TMP.sets/%f && rm $TMP.sets/%f || true"
done      
echo "" > /dev/stderr


rm $TMP*
