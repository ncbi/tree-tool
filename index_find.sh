#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 5 ]; then
  echo "Find approximately closest objects ordered by proximity"
  echo "An object is a list of attributes"
  echo "#1: directory with objects"
  echo "#2: #1 is large (0/1)"
  echo "#3: directory with attribute index"
  echo "#4: object size limit (S)"
  echo "#5: target object"  
  echo "Average time: O(A^2 S (A + log A + log N)), where N is the number of objects, and A is the average number of attributes in an object"
  exit 1
fi
OBJ_DIR=$1
LARGE=$2
INDEX=$3
SIZE_MAX=$4
TARGET=$5


if [ ! -e $TARGET ]; then
  error "$TARGET does not exist"
fi


TMP=`mktemp`
#echo $TMP
#set -x


# |$INDEX| + O(N A)

# $TMP.obj, $TMP.big_attr
# |$TMP.obj| = O(A S)
# Time = O(A (S + log N + log A))
touch $TMP.obj
while read ATTR
do
  if [ ! -e $INDEX/$ATTR ]; then  
    continue
  fi
  SIZE=`ls -l $INDEX/$ATTR | cut -f 5 -d ' '`
  if [ $SIZE -ge $SIZE_MAX ]; then
    echo $ATTR >> $TMP.big_attr
  else
    cat $INDEX/$ATTR >> $TMP.obj
  fi
done < $TARGET

# Time: O(A S (log A + log S))
sort -u $TMP.obj > $TMP.small_obj

# append: $TMP.obj
# Time: O(A^2 S (A + log N))
while read OBJ
do
  # Time: O(A (A + log N))
  H=""
  if [ $LARGE == 1 ]; then
    H=`$THIS/file2hash $OBJ`
  fi
  while read ATTR
  do
    if grep -x "$ATTR" $OBJ_DIR/$H/$OBJ &> /dev/null; then
      echo $OBJ >> $TMP.obj
    fi
  done < $TMP.big_attr
done < $TMP.small_obj
# |$TMP.obj| = O(A^2 S)

# Time: O(A^2 S (log A + log S))
sort $TMP.obj | uniq -c | sort -k 1 -n -r | sed 's/^ *[0-9]\+ \+//1'


rm $TMP*
