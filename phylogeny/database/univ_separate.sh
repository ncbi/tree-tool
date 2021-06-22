#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
GENOME="genome"
if [ $# -ne 3 ]; then
  echo "Create #1-univ-separate.dm"
  echo "Input: #1.list, $GENOME/*/[*/]*.prot-univ, hmm-univ.list"
  echo "Temporary output: #1.pairs.dir/, #1.pairs.dir.log/, #1-univ/"
  echo "#1: list of objects (file prefix)"
  echo "#2: 1 - BLOSUM62, 0 - PAM30"
  echo "#3: $GENOME/ is large (0/1)"
  exit 1
fi
F=$1
BLOSUM62=$2
LARGE=$3


LIST=$F.list
if [ ! -e $LIST ]; then
  error "File $LIST does not exist"
fi


TMP=`mktemp`
echo $TMP


# $F.pairs.dir
$THIS/../../list2pairs $LIST > $TMP
mkdir $F.pairs.dir
$THIS/../../splitList $TMP 500 $F.pairs.dir
  # PAR
  
mkdir $F-univ
$THIS/../../dir2log.sh $F.pairs.dir
$THIS/../../grid_wait.sh 1
$THIS/../../trav $F.pairs.dir.log "$QSUB_5 -N j%f '$THIS/prot_collection2dissim_separate.sh %f $BLOSUM62 $F.pairs.dir $F-univ $LARGE' > /dev/null"
$THIS/../../qstat_wait.sh 2000 1
rmdir $F.pairs.dir.log/

rm -r $F.pairs.dir &

$THIS/../../trav  $F-univ "cat %d/%f" > $TMP.dissim
rm -r $F-univ/ &

N=`cat $TMP.dissim | wc -l`
echo "OBJNUM $N name nomult"                >  $F-univ-separate.dm
echo "ATTRIBUTES"                           >> $F-univ-separate.dm
cat hmm-univ.list | sed 's/$/ Positive 6/1' >> $F-univ-separate.dm
echo "DATA"                                 >> $F-univ-separate.dm
cat $TMP.dissim | sed -e 's/\t/-/1'         >> $F-univ-separate.dm


rm $TMP*
wait
