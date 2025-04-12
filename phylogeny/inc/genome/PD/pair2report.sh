#!/bin/bash --noprofile
#source bash_common.sh
if [ $# -ne 3 ]; then
  echo "Print dissimilarity statistics: # identical proteins, # common universal proteins, AAI%"
  echo "#1: object 1 or new object path, where the object name is basename"
  echo "#2: #1 is new (0/1)"
  echo "#3: object 2"
  exit 1
fi
OBJ1=$1
NEW=$2
OBJ2=$3


error "$0 is not implemented"
