#!/bin/bash --noprofile
THIS=$( realpath "$( dirname $0 )" )
source $THIS/bash_common.sh


super_section "numeric"
$THIS/numeric_test -qc go 

super_section "graph"
$THIS/graph_test go -qc

super_section "assignment"
N=$( $THIS/assignment -qc $THIS/assignment.txt -noprogress )
M=$( $THIS/assignment -qc $THIS/assignment.txt -noprogress -verbose 1 | grep "total_lo" | tr '\t' ' ' | cut -f 3 -d ' ' | count | grep -w "^sum" | cut -f 2 )
if [ $N -lt $M ]; then
  error "assigment: $N < $M"
fi

super_section "dm"
# Time: 31 sec.
$THIS/dm/dm_test.sh 1

super_section "kmerIndex"
# Time: 18 sec.
$THIS/genetics/kmerIndex_test.sh go

super_section "distTree"
# Time: 75 (63) min.
time $THIS/phylogeny/distTree_test.sh go

super_section "featureTree"
# Time: 10 min.
time $THIS/phylogeny/featureTree_test.sh go

super_section "xml"
$THIS/xml/xml_test.sh 0


success
