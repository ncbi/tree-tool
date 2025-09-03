#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
echo "Compare the files of two directories, except .-files and binary files"
  echo "#1: file/directory 1"
  echo "#2: file/directory 2"
  exit 1
fi
D1=$1
D2=$2


TMP=$( mktemp )


if [ -d $D1 ] && [ ! -d $D2 ]; then
  warning "$D1 is a directory, but $D2 is not"
elif [ ! -d $D1 ] && [ -d $D2 ]; then
  warning "$D2 is a directory, but $D1 is not"
elif [ -d $D1 ] && [ -d $D2 ]; then
  ls $D1 > $TMP.1
  ls $D2 > $TMP.2

  $THIS/setMinus $TMP.1 $TMP.2 > $TMP.1-2
  if [ -s $TMP.1-2 ]; then
    warning "In $D1, but not in $D2:"
    cat $TMP.1-2
  fi

  $THIS/setMinus $TMP.2 $TMP.1 > $TMP.2-1
  if [ -s $TMP.2-1 ]; then
    warning "In $D2, but not in $D1:"
    cat $TMP.2-1
  fi

  $THIS/setIntersect.sh $TMP.1 $TMP.2 0 > $TMP.common
  $THIS/trav -step 1  -noprogress  $TMP.common "$THIS/dir_diff.sh $D1/%f $D2/%f"
else
  file -b $D1 > $TMP.1
  file -b $D2 > $TMP.2
  if grep -wq "ASCII" $TMP.1 && grep -wq "ASCII" $TMP.2; then
    diff $D1 $D2 > $TMP || true
    if [ -s $TMP ]; then
      warning "diff $D1 $D2"
    fi
  elif diff $TMP.1 $TMP.2 ; then
    :
  else
    warning "$D1 and $D2 have different file types"
  fi
fi


rm $TMP*
