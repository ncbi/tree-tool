#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Find the place of an ITS sequence (ITS1 - 5.8S rRNA - ITS2) in the Fungal ITS tree"
  echo "#1: ITS sequence"
  echo "#2: ITS directory"
  exit 1
fi
ACC=$1
DIR=$2  # /home/brovervv/panfs/marker/Fungi/ITS


TMP=`mktemp`  


NAME=`grep ">" $ACC | sed 's/^>//1' | sed 's/ .*$//1'`

blastn  -db $DIR/inc/seq.fa  -query $1  -show_gis  -evalue 1e-20  -outfmt '6 sseqid' > $TMP.blastn
head -100 $TMP.blastn | sed 's/$/ '$NAME'/1' > $TMP.request

cp /dev/null $TMP.dissim
echo "FAIL" > $TMP.leaf
while [ -s $TMP.request ]; do
  $THIS/../dissim/dna_pair2dissim  -coeff 0.0082  -name_new $NAME  -file_new $ACC  $TMP.request $DIR/seq 140 $TMP.dissim-add &> /dev/null
  rm $TMP.request
  cat $TMP.dissim-add >> $TMP.dissim
  $THIS/../phylogeny/distTree_new  $DIR/inc/  -variance lin  -name $NAME  -dissim $TMP.dissim  -request $TMP.request  -leaf $TMP.leaf
done

L=(`cat $TMP.leaf`)
echo "${L[0]} has arc of length ${L[2]} joining above ${L[1]} by ${L[3]} "


rm -fr $TMP*

