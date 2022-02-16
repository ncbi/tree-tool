#!/bin/bash --noprofile
THIS=`dirname $0`
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


TMP=`mktemp`
#echo $TMP > /dev/stderr 


# PAR
T=90
LEN_MIN=100


$THIS/dna_coverage.sh $DNA $PL combine $T $LEN_MIN '' 0 | tail -n +2 | awk '$6 >= '$T' && $7 >= '$T | cut -f 2,10 > $TMP.contig-pls

cut -f 1 $TMP.contig-pls > $TMP.contig
$THIS/extractFastaDna $DNA $TMP.contig -noprogress > $TMP.contig-fa

cut -f 2 $TMP.contig-pls | tr ' ' '\n' | sort -u > $TMP.pls
$THIS/../trav -noprogress $TMP.pls "$THIS/fasta2len -noprogress $PL/%f" > $TMP.len
sort -k2,2nr $TMP.len | cut -f 1 > $TMP.pls-sorted
PLS=(`cat $TMP.pls-sorted`)

echo -e "#Plasmid\tlen\tcoverage\tpcoverage\tpident\tContigs" > $TMP.res
i=0
while [ $i -lt ${#PLS[@]} ]; do
  $THIS/dna_coverage.sh $PL/${PLS[i]} $TMP.contig-fa combine $T $LEN_MIN '' 0 > $TMP.pl-cov
  tail -n +2 $TMP.pl-cov | awk '$6 >= '$T' && $7 >= '$T | cut -f 10 | tr ' ' '\n' | sort > $TMP.contigs
  if [ -s $TMP.contigs ]; then
    tail -n +2 $TMP.pl-cov | cut -f 1,4,5,6,7,10 >> $TMP.res
    $THIS/extractFastaDna $TMP.contig-fa $TMP.contigs -remove -noprogress > $TMP.contig-fa-new
    mv $TMP.contig-fa-new $TMP.contig-fa
  fi
  i=$(( $i + 1 ))
done

cat $TMP.res


rm $TMP*  
