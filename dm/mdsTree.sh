#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Hierarchical clustering by MDS. Output directory: #1/"
  echo "#1: Input .dm file without .dm"
  echo "#2: Pairwise measurement attribute in #1"
  echo "#3: mds -attrType"
  exit 1
fi


DEBUG=0  # or 1 


TMP=`mktemp` 


# TMP, VERBOSE
if [ $DEBUG == 1 ]; then
  VERBOSE=2
  echo $TMP
  set -x
else
  VERBOSE=0
fi

EXITCODE=1


mkdir $1.dir

rm -rf $TMP.tmp
mkdir $TMP.tmp

cp $1.dm $TMP.tmp/

echo "" > $TMP.mds

while true; do
  echo ""
  F_FULL=`ls $TMP.tmp`
  if [ -z "$F_FULL" ]; then
    break
  fi
  F=($F_FULL)
  F=`basename ${F[0]} .dm`
  grep -v '^#' $TMP.tmp/$F.dm > $TMP.grep || true
  N=(`head -1 $TMP.grep`)
  set -o errexit
  echo "# objects =  ${N[1]}"
  if [ ${N[1]} -le 5 ]; then  # PAR
    mv $TMP.tmp/$F.dm $1.dir/
  else
    echo "" >> $TMP.mds
    echo "" >> $TMP.mds
    L_old=`ls $TMP.tmp/ | wc -l`
    $THIS/mds -verbose $VERBOSE  -attrType $3  -maxClusters 6  -clusteringDir $TMP.tmp  -attr $2  $TMP.tmp/$F >> $TMP.mds
    L_new=`ls $TMP.tmp/ | wc -l`
    if [ $L_old == $L_new ]; then
      mv $TMP.tmp/$F.dm $1.dir/
    else
      rm $TMP.tmp/$F.dm
    fi
  fi
done


if [ $DEBUG == 0 ]; then
  rm -fr $TMP*   
fi
