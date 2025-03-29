#!/bin/bash --noprofile
THIS=$( dirname $0 )
source CPP_DIR/bash_common.sh
if [ $# -ne 2 ]; then
  echo "QC an object"
  echo "#1: object prefix"
  echo "#2: verbose: 0/1"
  exit 1
fi
PREF=$1
VERB=$2


if [ $VERB == 1 ]; then
  set -x
fi


function check 
{
  local SUF=$1
  local FILE="$2"
  local F=$PREF.$SUF
  CPP_DIR/check_file.sh $F 1
  if [ "$FILE" == "gzip compressed data" ]; then
    file $F | grep -q ": $FILE"
  else
    if [ -s $F ]; then
      file -m /dev/null $F | grep -q ": $FILE"
    fi
  fi
}

#check "hash-CDS"  ": ASCII text" || check "hash-CDS"  ": symbolic link to "  # File is not used
check "hash-HMM"       "ASCII text"
check "hash-PRT"       "ASCII text"
check "univ"           "ASCII text"
check "prot-univ"      "symbolic link to "
check "prot-univ.HMM"  "ASCII text"
check "prot.gz"        "gzip compressed data"
check "stat"           "ASCII text"

if [ -s $PREF.prot-univ.HMM ]; then
  head -1 $PREF.prot-univ.HMM | grep -vq " ref="
fi

if [ -s $PREF.prot-univ.tblastn ]; then
  head -1 $PREF.prot-univ.tblastn | grep -q " ref="
fi

PU=$PREF.prot-univ
REF=$( basename $( realpath $PU ) )
if [ $REF != $( basename $PU.HMM ) ]; then
  # && [ $REF != $( basename $PU.tblastn ) ]
  error "Strange reference of $PU: $REF"
fi
