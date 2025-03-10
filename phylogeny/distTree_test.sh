#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: 70 min."
  exit 1
fi


TMP=$( mktemp )
if [ ${AT_NCBI:-0} == 1 ]; then
  comment $TMP
fi


DATA=$THIS/data
#set -x


#if false; then  
section "mdsTree: Enterobacteriaceae"
#rm -rf $DATA/Enterobacteriaceae.dir/
gunzip -c $DATA/Enterobacteriaceae.dm.gz > $TMP.dm
$THIS/../dm/mdsTree.sh $TMP "Conservation" 2  &> $TMP.out
$THIS/makeDistTree  -qc -input_tree $TMP.dir/  -data $TMP  -dissim_attr "Conservation"  -variance linExp  -optimize  -output_tree $TMP.Enterobacteriaceae.tree > $TMP.Enterobacteriaceae.distTree
$THIS/distTree_compare_criteria.sh $TMP.Enterobacteriaceae.distTree $DATA/Enterobacteriaceae.distTree

if false; then
  section "To Newick"
  gunzip -c $DATA/Enterobacteriaceae.dm.gz > $TMP.dm
  $THIS/printDistTree  -qc  -data $TMP  -dissim_attr "Conservation"  -variance linExp  $TMP.Enterobacteriaceae.tree  -order  -decimals 4  -ext_name > $TMP.Enterobacteriaceae.nw
  diff $TMP.Enterobacteriaceae.nw $DATA/Enterobacteriaceae.nw
fi

section "From Newick"
$THIS/newick2tree -qc $DATA/Enterobacteriaceae.nw > $TMP.Enterobacteriaceae.tree
$THIS/printDistTree -qc $TMP.Enterobacteriaceae.tree  -decimals 4  -ext_name > $TMP.Enterobacteriaceae.nw
diff $TMP.Enterobacteriaceae.nw $DATA/Enterobacteriaceae.nw

section "Perfect tree"
$THIS/makeDistTree  -qc  -data $DATA/tree4  -variance linExp  -optimize  -output_tree $TMP.tree4 | grep -v '^CHRON: ' > $TMP.tree4.makeDistTree
$THIS/distTree_compare_criteria.sh $TMP.tree4.makeDistTree $DATA/tree4.makeDistTree
$THIS/statDistTree $TMP.tree4 > $TMP.tree4.stat
diff $TMP.tree4.stat $DATA/tree4.stat
comment "Verbose"
$THIS/makeDistTree -qc  -data $DATA/tree4  -variance linExp  -optimize  -verbose 2 &> $TMP.out

section "Random tree"
gunzip -c $DATA/randomTree.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance lin  -output_tree $TMP.random-output.tree > $TMP.out
$THIS/makeDistTree  -qc  -input_tree $TMP.random-output.tree    -data $TMP  -variance lin  -optimize | grep -v '^CHRON: ' > $TMP.randomTree.makeDistTree
$THIS/distTree_compare_criteria.sh $TMP.randomTree.makeDistTree $DATA/randomTree.makeDistTree
comment "Verbose"
$THIS/makeDistTree  -qc  -input_tree $TMP.random-output.tree    -data $TMP  -variance lin  -verbose 2 &> $TMP.out

section "Salmonella"
# Check time ??
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 10 \
  -delete_criterion_outliers $TMP.Salmonella.criterion_outliers \
  -delete_deformation_outliers $TMP.Salmonella.deformation_outliers \
  -delete_hybrids $TMP.Salmonella.hybrids \
  | grep -v '^CHRON: ' > $TMP.Salmonella.distTree
diff $TMP.Salmonella.criterion_outliers $DATA/Salmonella.criterion_outliers
$THIS/../sort.sh $TMP.Salmonella.deformation_outliers
cut -f 1  $TMP.Salmonella.deformation_outliers > $TMP.Salmonella.deformation_outliers1
cut -f 1 $DATA/Salmonella.deformation_outliers > $TMP.data.Salmonella.deformation_outliers1
diff $TMP.Salmonella.deformation_outliers1 $TMP.data.Salmonella.deformation_outliers1
diff $TMP.Salmonella.hybrids $DATA/Salmonella.hybrids
$THIS/distTree_compare_criteria.sh $TMP.Salmonella.distTree $DATA/Salmonella.distTree

section "testDistTree"
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -data $TMP  -variance linExp  -optimize  -output_tree $TMP.tree > $TMP.out
$THIS/testDistTree -qc  $TMP  -input_tree  $TMP.tree  -variance linExp 

section "Salmonella: delete"
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -delete $DATA/delete.list  -check_delete > $TMP.out

if [ -d $DATA/inc.ITS ]; then
  section "ITS threads"
  $THIS/makeDistTree  -threads 10  -data $DATA/inc.ITS/  -variance linExp  -variance_dissim  -optimize  -skip_len  -reinsert  -subgraph_iter_max 1  -noqual > $TMP.ITS.distTree
  $THIS/distTree_compare_criteria.sh $TMP.ITS.distTree $DATA/ITS.distTree
  
  section "distTree_new"
  $THIS/distTree_new -qc $DATA/inc.ITS/  -variance linExp
  diff $DATA/inc.ITS/search/NR_073289.1/leaf    $DATA/inc.ITS/leaf.expected
  diff $DATA/inc.ITS/search/NR_073289.1/request $DATA/inc.ITS/request.expected
  rm $DATA/inc.ITS/search/NR_073289.1/leaf
  rm $DATA/inc.ITS/search/NR_073289.1/request
fi

section "Saccharomyces hybrids"
gunzip -c $DATA/Saccharomyces.dm.gz > $TMP.dm
$THIS/makeDistTree -qc  -threads 10  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 2  \
  -hybridness_min 1.2  -hybrid_parent_pairs $TMP.Saccharomyces.hybrid_parent_pairs  -delete_hybrids $TMP.Saccharomyces.hybrid  -dissim_boundary 0.675 \
  -delete_criterion_outliers $TMP.Saccharomyces.criterion_outliers  -criterion_outlier_num_max 1 \
  -delete_deformation_outliers $TMP.Saccharomyces.deformation_outliers  -deformation_outlier_num_max 1 \
  -output_tree $TMP.tree > $TMP.out
diff $TMP.Saccharomyces.hybrid $DATA/Saccharomyces.hybrid
$THIS/hybrid2list.sh $TMP.Saccharomyces.hybrid > $TMP.Saccharomyces.hybrid.list
diff $TMP.Saccharomyces.hybrid_parent_pairs $DATA/Saccharomyces.hybrid_parent_pairs
diff $TMP.Saccharomyces.criterion_outliers $DATA/Saccharomyces.criterion_outliers
cut -f 1  $TMP.Saccharomyces.deformation_outliers > $TMP.Saccharomyces.deformation_outliers1
cut -f 1 $DATA/Saccharomyces.deformation_outliers > $TMP.data.Saccharomyces.deformation_outliers1
diff $TMP.Saccharomyces.deformation_outliers1 $TMP.data.Saccharomyces.deformation_outliers1
# Saccharomyces.distTree
$THIS/tree2obj.sh $TMP.tree > $TMP.list
$THIS/../dm/dm2subset $TMP $TMP.list > $TMP.subset.dm
$THIS/makeDistTree  -threads 10  -data $TMP.subset  -input_tree $TMP.tree  -variance linExp  -optimize  -reinsert  -subgraph_iter_max 10  > $TMP.Saccharomyces.distTree
$THIS/distTree_compare_criteria.sh $TMP.Saccharomyces.distTree $DATA/Saccharomyces.distTree

section "-variance_min"
# 0.0005 = average arc length / 100
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -variance_min 0.0005  -optimize  -subgraph_iter_max 10 \
  -delete_criterion_outliers $TMP.Salmonella-var_min.criterion_outliers \
  -delete_deformation_outliers $TMP.Salmonella-var_min.deformation_outliers \
  -delete_hybrids $TMP.Salmonella-var_min.hybrids \
  | grep -v '^CHRON: ' > $TMP.Salmonella-var_min.distTree
diff $TMP.Salmonella-var_min.criterion_outliers $DATA/Salmonella-var_min.criterion_outliers
$THIS/../sort.sh $TMP.Salmonella-var_min.deformation_outliers 
cut -f 1  $TMP.Salmonella-var_min.deformation_outliers > $TMP.Salmonella-var_min.deformation_outliers1
cut -f 1 $DATA/Salmonella-var_min.deformation_outliers > $TMP.data.Salmonella-var_min.deformation_outliers1
diff $TMP.Salmonella-var_min.deformation_outliers1 $TMP.data.Salmonella-var_min.deformation_outliers1
diff $TMP.Salmonella-var_min.hybrids $DATA/Salmonella-var_min.hybrids
$THIS/distTree_compare_criteria.sh $TMP.Salmonella-var_min.distTree $DATA/Salmonella-var_min.distTree


if [ ${AT_NCBI:-0} == 0 ]; then
  rm Saccharomyces.hybrid.list
else
  echo ""
  super_section "Two dissimilarity types"

  section "Salmonella2"
  gunzip -c $DATA/Salmonella2.dm.gz > $TMP.dm
  $THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 10  \
    -output_dissim_coeff $TMP.Salmonella2.dissim_coeff \
    -delete_criterion_outliers $TMP.Salmonella2.criterion_outliers \
    -delete_deformation_outliers $TMP.Salmonella2.deformation_outliers \
    -delete_hybrids $TMP.Salmonella2.hybrids \
    | grep -v '^CHRON: ' > $TMP.Salmonella2.distTree
  diff $TMP.Salmonella2.dissim_coeff $DATA/Salmonella2.dissim_coeff
  diff $TMP.Salmonella2.criterion_outliers $DATA/Salmonella2.criterion_outliers
  $THIS/../sort.sh $TMP.Salmonella2.deformation_outliers
  cut -f 1  $TMP.Salmonella2.deformation_outliers > $TMP.Salmonella2.deformation_outliers1
  cut -f 1 $DATA/Salmonella2.deformation_outliers > $TMP.data.Salmonella2.deformation_outliers1
  diff $TMP.Salmonella2.deformation_outliers1 $TMP.data.Salmonella2.deformation_outliers1
  diff $TMP.Salmonella2.hybrids $DATA/Salmonella2.hybrids
  $THIS/distTree_compare_criteria.sh $TMP.Salmonella2.distTree $DATA/Salmonella2.distTree

  section "Saccharomyces hybrids"
  gunzip -c $DATA/Saccharomyces2.dm.gz > $TMP.dm
  $THIS/makeDistTree -qc  -threads 10  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 2  \
    -hybridness_min 1.2  -hybrid_parent_pairs $TMP.Saccharomyces2.hybrid_parent_pairs  -delete_hybrids $TMP.Saccharomyces2.hybrid  -dissim_boundary 0.675 \
    -delete_criterion_outliers $TMP.Saccharomyces2.criterion_outliers  -criterion_outlier_num_max 1 \
    -delete_deformation_outliers $TMP.Saccharomyces2.deformation_outliers  -deformation_outlier_num_max 1 \
    -output_tree $TMP.tree > $TMP.out
  $THIS/../tsv/tsv_cut.sh $TMP.Saccharomyces2.hybrid "1,2,3,4,5,6,7,8,9" 0 | sort > $TMP.Saccharomyces2.hybrid.cut
  diff $TMP.Saccharomyces2.hybrid.cut $DATA/Saccharomyces2.hybrid
  $THIS/hybrid2list.sh $TMP.Saccharomyces2.hybrid > $TMP.Saccharomyces2.hybrid.list
  $THIS/../tsv/tsv_cut.sh $TMP.Saccharomyces2.hybrid_parent_pairs "1,2,3,4,5,6,7,8,9,10,11,12" 0 | sort > $TMP.Saccharomyces2.hybrid_parent_pairs.cut
  diff -b $TMP.Saccharomyces2.hybrid_parent_pairs.cut $DATA/Saccharomyces2.hybrid_parent_pairs
  diff $TMP.Saccharomyces2.criterion_outliers $DATA/Saccharomyces2.criterion_outliers
  cut -f 1  $TMP.Saccharomyces2.deformation_outliers > $TMP.Saccharomyces2.deformation_outliers1
  cut -f 1 $DATA/Saccharomyces2.deformation_outliers > $TMP.data.Saccharomyces2.deformation_outliers1
  diff $TMP.Saccharomyces2.deformation_outliers1 $TMP.data.Saccharomyces2.deformation_outliers1
  # Saccharomyces2.distTree
  $THIS/tree2obj.sh $TMP.tree > $TMP.list
  $THIS/../dm/dm2subset $TMP $TMP.list > $TMP.subset.dm
  $THIS/makeDistTree  -threads 10  -data $TMP.subset  -input_tree $TMP.tree  -variance linExp  -optimize  -reinsert  -subgraph_iter_max 10  > $TMP.Saccharomyces2.distTree
  $THIS/distTree_compare_criteria.sh $TMP.Saccharomyces2.distTree $DATA/Saccharomyces2.distTree

  set +o errexit
  N=$( diff -y --suppress-common-lines Saccharomyces.hybrid.list Saccharomyces2.hybrid.list | wc -l )
  set -o errexit
  if [ $N -gt 1 ]; then  # ??
    diff $TMP.Saccharomyces.hybrid.list $TMP.Saccharomyces2.hybrid.list
  fi
  diff $DATA/Saccharomyces2.criterion_outliers $DATA/Saccharomyces.criterion_outliers
  cut -f 1 $DATA/Saccharomyces2.deformation_outliers > $TMP.Saccharomyces2.deformation_outliers1
  cut -f 1 $DATA/Saccharomyces.deformation_outliers  > $TMP.Saccharomyces.deformation_outliers1
  diff $TMP.Saccharomyces2.deformation_outliers1 $TMP.Saccharomyces.deformation_outliers1


  super_section "Many dissimilarity types"

  section "Virus9"
  gunzip -c $DATA/Virus9.dm.gz > $TMP.dm
  $THIS/makeDistTree  -qc  -data $TMP  -optimize  -subgraph_iter_max 100  -variance linExp  -variance_dissim  -output_dissim_coeff $TMP.Virus9.coeff  -output_data $TMP.Virus9-out  -output_tree $TMP.Virus9.tree  1> $TMP.1 2> $TMP.out
  diff $TMP.Virus9.coeff $DATA/Virus9.coeff
  A=$( grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1 )
  #
  $THIS/makeDistTree  -qc  -data $TMP.Virus9-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree $TMP.Virus9-out.tree  1> $TMP.2 2> $TMP.out
  B=$( grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1 )
  if [ "$A" != "$B" ]; then
    error "$A != $B"
  fi
  $THIS/printDistTree  -qc  $TMP.Virus9.tree      -order  -decimals 3  > $TMP.Virus9.nw
  $THIS/printDistTree  -qc  $TMP.Virus9-out.tree  -order  -decimals 3  > $TMP.Virus9-out.nw
  diff $TMP.Virus9.nw $TMP.Virus9-out.nw
  

  section "Virus110"
  gunzip -c $DATA/Virus110.dm.gz > $TMP.dm
  $THIS/makeDistTree  -data $TMP  -variance_min 0.005  -variance linExp  -variance_dissim  -optimize  -subgraph_iter_max 100  -output_data $TMP.Virus110-out  -output_tree $TMP.Virus110.tree 1> $TMP.1 
  A=$( grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1 )
  #
  $THIS/makeDistTree  -qc  -input_tree $TMP.Virus110.tree  -data $TMP.Virus110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree $TMP.Virus110-out1.tree 1> $TMP.2 2> $TMP.out
  B=$( grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1 )
  if [ "$A" != "$B" ]; then
    error "$A != $B"
  fi
  $THIS/printDistTree  -qc  $TMP.Virus110.tree       -order  -decimals 1 | sed 's/,(/,\n(/g' > $TMP.Virus110.nw
  $THIS/printDistTree  -qc  $TMP.Virus110-out1.tree  -order  -decimals 1 | sed 's/,(/,\n(/g' > $TMP.Virus110-out1.nw
  diff $TMP.Virus110.nw $TMP.Virus110-out1.nw || true
  #
  if false; then  # ??
    makeDistTree  -qc  -data $TMP.Virus110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree $TMP.Virus110-out2.tree 1> $TMP.2 2> $TMP.out
    B=$( grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1 )
    if [ "$A" != "$B" ]; then
      error "$A != $B"
    fi
    $THIS/printDistTree  -qc  $TMP.Virus110.tree       -order  -decimals 2  > $TMP.Virus110.nw
    $THIS/printDistTree  -qc  $TMP.Virus110-out2.tree  -order  -decimals 2  > $TMP.Virus110-out2.nw
    diff $TMP.Virus110.nw $TMP.Virus110-out.nw
  fi
  #
fi


rm -r $TMP*


success
