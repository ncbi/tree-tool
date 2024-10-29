#!/bin/bash --noprofile
THIS=$( realpath "$( dirname $0 )" )
source $THIS/bash_common.sh


super_section "graph"
$THIS/graph_test go -qc

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
