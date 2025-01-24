#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# != 2 ]; then
  echo "Plasmid finder"
  echo "Print plasmid report"
  echo "#1: DNA FASTA file"
  echo "#2: BLAST database of plasmids and directory with separate plasmids"
  exit 1
fi
DNA=$1
PL=$2


TMP=$( mktemp )
#comment $TMP 


# PAR
IDENT=98
LEN=100
CONTIG_COV=90
PLASMID_COV=98


$THIS/dna_coverage.sh $DNA $PL combine $IDENT $LEN '' 0 0 | tail -n +2 | awk '$5 >= '$CONTIG_COV | cut -f 2,14 > $TMP.contig-pls

cut -f 1 $TMP.contig-pls > $TMP.contig
#$THIS/extractFastaDna $DNA $TMP.contig -noprogress > $TMP.contig-fa
$THIS/filterFasta $DNA  -target $TMP.contig  -noprogress > $TMP.contig-fa

cut -f 2 $TMP.contig-pls | tr ' ' '\n' | sort -u > $TMP.pls
$THIS/../trav -noprogress $TMP.pls "$THIS/fasta2len -noprogress $PL/%f" > $TMP.len
sort -k2,2nr $TMP.len | cut -f 1 > $TMP.pls-sorted
PLS=( $( cat $TMP.pls-sorted ) )

echo -e "#plasmid\tcontig" 
i=0
while [ $i -lt ${#PLS[@]} ]; do
  $THIS/dna_coverage.sh $PL/${PLS[i]} $TMP.contig-fa combine $IDENT $LEN '' 0 0 | tail -n +2 | awk '$5 >= '$PLASMID_COV > $TMP.pl-contigs-$i
  cut -f 14 $TMP.pl-contigs-$i | tr ' ' '\n' | sort > $TMP.contigs-$i
  if [ -s $TMP.contigs-$i ]; then
   #$THIS/extractFastaDna $TMP.contig-fa $TMP.contigs-$i -noprogress > $TMP.pl-contig.fa-$i
    $THIS/filterFasta $TMP.contig-fa  -target $TMP.contigs-$i  -noprogress > $TMP.pl-contig.fa-$i
    $THIS/dna_coverage.sh $TMP.pl-contig.fa-$i $PL/${PLS[i]} combine $IDENT $LEN '' 0 0 | tail -n +2 | awk '$5 >= '$CONTIG_COV | cut -f 2 > $TMP.contigs-match-$i
    sed 's/^/'${PLS[i]}'\t/1' $TMP.contigs-match-$i 
   #$THIS/extractFastaDna $TMP.contig-fa $TMP.contigs-match-$i -remove -noprogress > $TMP.contig-fa-new-$i
    $THIS/filterFasta $TMP.contig-fa  -target $TMP.contigs-match-$i  -remove -noprogress > $TMP.contig-fa-new-$i
    mv $TMP.contig-fa-new-$i $TMP.contig-fa
  fi
  i=$(( i + 1 ))
done


rm $TMP*  

