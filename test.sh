#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/bash_common.sh


super_section "graph"
$THIS/graph_test go -qc

super_section "dm"
$THIS/dm/dm_test.sh 1

super_section "kmerIndex"
$THIS/genetics/kmerIndex_test.sh go

super_section "distTree"
$THIS/phylogeny/distTree_test.sh go

super_section "featureTree"
$THIS/phylogeny/featureTree_test.sh go

super_section "xml"
$THIS/xml/xml_test.sh 0
