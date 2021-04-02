#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 4 ]; then
  echo "Find approximately closest objects ordered by proximity"
  echo "An object is a list of attributes"
  echo "#1: directory with objects"
  echo "#2: directory with attribute index"
  echo "#3: object size limit (S)"
  echo "#4: target object"  
  echo "Average time: O(S A^2 (A + log S)), where A is the average number of attributes in an object"
  exit 1
fi
OBJ_DIR=$1
INDEX=$2
SIZE_MAX=$3
TARGET=$4


TMP=`mktemp`
#echo $TMP


# $TMP.obj, $TMP.big_attr
# |$TMP.obj| = O(S A)
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

# Time: O(S A (log A + log A))
sort -u $TMP.obj > $TMP.small_obj

# append: $TMP.obj
# Time: O(S A^3)
while read OBJ
do
  # Time: O(A^2)
  while read ATTR
  do
    if grep -x "$ATTR" $OBJ_DIR/$OBJ &> /dev/null; then
      echo $OBJ >> $TMP.obj
    fi
  done < $TMP.big_attr
done < $TMP.small_obj
# |$TMP.obj| = O(S A^2)

# Time: O(S A^2 (log S + log A))
sort $TMP.obj | uniq -c | sort -k 1 -n -r | sed 's/^ *[0-9]\+ \+//1'


rm $TMP*
