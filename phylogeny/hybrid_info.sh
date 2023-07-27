#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Information on a hybrid object"
  echo "#1: incremental distance tree directory"
  echo "#2: object id or external identifier"
  exit 1
fi
INC=$1
NAME=$2


META=$INC/../Metadata.tsv
if [ ! -e $META ]; then
  error "Metadata file is not found"
fi


TMP=`mktemp`
#comment $TMP 


tail -n +2 $META | cut -f 1,2 > $TMP.dict


# OBJ
grep -w $NAME $TMP.dict > $TMP || true
if [ ! -s $TMP ]; then
  error "Object $NAME is not found in $META"
fi
N=`cat $TMP | wc -l`
if [ $N -gt 1 ]; then
  cat $TMP
  error "Object $NAME is not unique"
fi
OBJ=`cut -f 1 $TMP`

grep -w $OBJ $INC/hist/hybrid.* > $TMP.grep || true
if [ ! -s $TMP.grep ]; then
  error "No hybrid information"
fi

echo -e "#iteration\tchild\thybridness\tparent1\tparent2\tdist_parent1\tdist_parent2" >  $TMP.tsv
sed 's/^[^.]*\.\([^:]\+\):/\1\t/1' $TMP.grep | sed 's/$/\t/1' | cut -f 1-7            >> $TMP.tsv
  #                                            ^ for an old format

echo -e "#child\tchild_acc" >  $TMP.child
cat $TMP.dict               >> $TMP.child
$THIS/../tsv/tsv_expand.sh $TMP.tsv $TMP.child '' &> /dev/null

echo -e "#parent1\tparent1_acc" >  $TMP.1
cat $TMP.dict                   >> $TMP.1
$THIS/../tsv/tsv_expand.sh $TMP.tsv $TMP.1 '' &> /dev/null

echo -e "#parent2\tparent2_acc" >  $TMP.2
cat $TMP.dict                   >> $TMP.2
$THIS/../tsv/tsv_expand.sh $TMP.tsv $TMP.2 '' &> /dev/null


cat $TMP.tsv


rm -rf $TMP*
