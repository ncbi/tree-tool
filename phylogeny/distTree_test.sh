#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: 120 min."
  exit 1
fi


TMP=`mktemp`
if [ ${AT_NCBI:-0} == 1 ]; then
  comment $TMP
fi


DATA=$THIS/data


#if false; then  
section "mdsTree: Enterobacteriaceae"
#rm -rf $DATA/Enterobacteriaceae.dir/
gunzip -c $DATA/Enterobacteriaceae.dm.gz > $TMP.dm
$THIS/../dm/mdsTree.sh $TMP Conservation 2  &> $TMP.out
$THIS/makeDistTree  -qc -input_tree $TMP.dir/  -data $TMP  -dissim_attr "Conservation"  -variance linExp  -optimize  -output_tree Enterobacteriaceae.tree > Enterobacteriaceae.distTree
$THIS/distTree_compare_criteria.sh Enterobacteriaceae.distTree $DATA/Enterobacteriaceae.distTree
rm Enterobacteriaceae.distTree
#rm -r $DATA/Enterobacteriaceae.dir/

section "To Newick"
gunzip -c $DATA/Enterobacteriaceae.dm.gz > $TMP.dm
$THIS/printDistTree  -qc  -data $TMP  -dissim_attr "Conservation"  -variance linExp  Enterobacteriaceae.tree  -order  -decimals 4  -ext_name > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw $DATA/Enterobacteriaceae.nw
rm Enterobacteriaceae.nw

rm Enterobacteriaceae.tree

section "From Newick"
$THIS/newick2tree -qc $DATA/Enterobacteriaceae.nw > Enterobacteriaceae.tree
$THIS/printDistTree -qc Enterobacteriaceae.tree  -decimals 4  -ext_name > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw $DATA/Enterobacteriaceae.nw
rm Enterobacteriaceae.nw
rm Enterobacteriaceae.tree

section "Perfect tree"
$THIS/makeDistTree  -qc  -data $DATA/tree4  -variance linExp  -optimize  -output_tree tree4 | grep -v '^CHRON: ' > tree4.makeDistTree
$THIS/distTree_compare_criteria.sh tree4.makeDistTree $DATA/tree4.makeDistTree
rm tree4.makeDistTree
$THIS/statDistTree tree4 > tree4.stat
rm tree4
diff tree4.stat $DATA/tree4.stat
rm tree4.stat
echo "Verbose ..."
$THIS/makeDistTree -qc  -data $DATA/tree4  -variance linExp  -optimize  -verbose 2 &> $TMP.out

section "Random tree"
gunzip -c $DATA/randomTree.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance lin  -output_tree random-output.tree > $TMP.out
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data $TMP  -variance lin  -optimize | grep -v '^CHRON: ' > randomTree.makeDistTree
$THIS/distTree_compare_criteria.sh randomTree.makeDistTree $DATA/randomTree.makeDistTree
echo "Verbose ..."
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data $TMP  -variance lin  -verbose 2 &> $TMP.out
rm randomTree.makeDistTree
rm random-output.tree

section "Salmonella"
# Check time ??
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 10 \
  -delete_criterion_outliers Salmonella.criterion_outliers \
  -delete_deformation_outliers Salmonella.deformation_outliers \
  -delete_hybrids Salmonella.hybrids \
  | grep -v '^CHRON: ' > Salmonella.distTree
diff Salmonella.criterion_outliers $DATA/Salmonella.criterion_outliers
rm Salmonella.criterion_outliers
$THIS/../sort.sh Salmonella.deformation_outliers
cut -f 1      Salmonella.deformation_outliers > $TMP.Salmonella.deformation_outliers
cut -f 1 $DATA/Salmonella.deformation_outliers > $TMP.data.Salmonella.deformation_outliers
diff $TMP.Salmonella.deformation_outliers $TMP.data.Salmonella.deformation_outliers
rm Salmonella.deformation_outliers
diff Salmonella.hybrids $DATA/Salmonella.hybrids
rm Salmonella.hybrids
$THIS/distTree_compare_criteria.sh Salmonella.distTree $DATA/Salmonella.distTree
rm Salmonella.distTree

section "testDistTree"
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -data $TMP  -variance linExp  -optimize  -output_tree $TMP.tree > $TMP.out
$THIS/testDistTree -qc  $TMP  -input_tree  $TMP.tree  -variance linExp 

section "Salmonella: delete"
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -delete $DATA/delete.list  -check_delete > $TMP.out

if [ -d $DATA/inc.ITS ]; then
  section "ITS threads"
  $THIS/makeDistTree  -threads 10  -data $DATA/inc.ITS/  -variance linExp  -variance_dissim  -optimize  -skip_len  -reinsert  -subgraph_iter_max 1  -noqual > ITS.distTree
  $THIS/distTree_compare_criteria.sh ITS.distTree $DATA/ITS.distTree
  rm ITS.distTree
fi

section "Saccharomyces hybrids"
gunzip -c $DATA/Saccharomyces.dm.gz > $TMP.dm
$THIS/makeDistTree -qc  -threads 10  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 2  \
  -hybridness_min 1.2  -hybrid_parent_pairs Saccharomyces.hybrid_parent_pairs  -delete_hybrids Saccharomyces.hybrid  -dissim_boundary 0.675 \
  -delete_criterion_outliers Saccharomyces.criterion_outliers  -criterion_outlier_num_max 1 \
  -delete_deformation_outliers Saccharomyces.deformation_outliers  -deformation_outlier_num_max 1 \
  -output_tree $TMP.tree > $TMP.out
diff Saccharomyces.hybrid $DATA/Saccharomyces.hybrid
$THIS/hybrid2list.sh Saccharomyces.hybrid > Saccharomyces.hybrid.list
rm Saccharomyces.hybrid
diff Saccharomyces.hybrid_parent_pairs $DATA/Saccharomyces.hybrid_parent_pairs
rm Saccharomyces.hybrid_parent_pairs
diff Saccharomyces.criterion_outliers $DATA/Saccharomyces.criterion_outliers
rm Saccharomyces.criterion_outliers
cut -f 1      Saccharomyces.deformation_outliers > $TMP.Saccharomyces.deformation_outliers
cut -f 1 $DATA/Saccharomyces.deformation_outliers > $TMP.data.Saccharomyces.deformation_outliers
diff $TMP.Saccharomyces.deformation_outliers $TMP.data.Saccharomyces.deformation_outliers
rm Saccharomyces.deformation_outliers
# Saccharomyces.distTree
$THIS/tree2obj.sh $TMP.tree > $TMP.list
$THIS/../dm/dm2subset $TMP $TMP.list > $TMP.subset.dm
$THIS/makeDistTree  -threads 10  -data $TMP.subset  -input_tree $TMP.tree  -variance linExp  -optimize  -reinsert  -subgraph_iter_max 10  > Saccharomyces.distTree
$THIS/distTree_compare_criteria.sh Saccharomyces.distTree $DATA/Saccharomyces.distTree
rm Saccharomyces.distTree

section "-variance_min"
# 0.0005 = average arc length / 100
gunzip -c $DATA/Salmonella.dm.gz > $TMP.dm
$THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -variance_min 0.0005  -optimize  -subgraph_iter_max 10 \
  -delete_criterion_outliers Salmonella-var_min.criterion_outliers \
  -delete_deformation_outliers Salmonella-var_min.deformation_outliers \
  -delete_hybrids Salmonella-var_min.hybrids \
  | grep -v '^CHRON: ' > Salmonella-var_min.distTree
diff Salmonella-var_min.criterion_outliers $DATA/Salmonella-var_min.criterion_outliers
rm Salmonella-var_min.criterion_outliers
$THIS/../sort.sh Salmonella-var_min.deformation_outliers 
cut -f 1      Salmonella-var_min.deformation_outliers > $TMP.Salmonella-var_min.deformation_outliers
cut -f 1 $DATA/Salmonella-var_min.deformation_outliers > $TMP.data.Salmonella-var_min.deformation_outliers
diff $TMP.Salmonella-var_min.deformation_outliers $TMP.data.Salmonella-var_min.deformation_outliers
rm Salmonella-var_min.deformation_outliers
diff Salmonella-var_min.hybrids $DATA/Salmonella-var_min.hybrids
rm Salmonella-var_min.hybrids
$THIS/distTree_compare_criteria.sh Salmonella-var_min.distTree $DATA/Salmonella-var_min.distTree
rm Salmonella-var_min.distTree


if [ ${AT_NCBI:-0} == 0 ]; then
  rm Saccharomyces.hybrid.list
else
  echo ""
  super_section "Two dissimilarity types"

  section "Salmonella2"
  gunzip -c $DATA/Salmonella2.dm.gz > $TMP.dm
  $THIS/makeDistTree  -qc  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 10  \
    -output_dissim_coeff Salmonella2.dissim_coeff \
    -delete_criterion_outliers Salmonella2.criterion_outliers \
    -delete_deformation_outliers Salmonella2.deformation_outliers \
    -delete_hybrids Salmonella2.hybrids \
    | grep -v '^CHRON: ' > Salmonella2.distTree
  diff Salmonella2.dissim_coeff $DATA/Salmonella2.dissim_coeff
  rm Salmonella2.dissim_coeff
  diff Salmonella2.criterion_outliers $DATA/Salmonella2.criterion_outliers
  rm Salmonella2.criterion_outliers
  $THIS/../sort.sh Salmonella2.deformation_outliers
  cut -f 1      Salmonella2.deformation_outliers > $TMP.Salmonella2.deformation_outliers
  cut -f 1 $DATA/Salmonella2.deformation_outliers > $TMP.data.Salmonella2.deformation_outliers
  diff $TMP.Salmonella2.deformation_outliers $TMP.data.Salmonella2.deformation_outliers
  rm Salmonella2.deformation_outliers
  diff Salmonella2.hybrids $DATA/Salmonella2.hybrids
  rm Salmonella2.hybrids
  $THIS/distTree_compare_criteria.sh Salmonella2.distTree $DATA/Salmonella2.distTree
  rm Salmonella2.distTree

  section "Saccharomyces hybrids"
  gunzip -c $DATA/Saccharomyces2.dm.gz > $TMP.dm
  $THIS/makeDistTree -qc  -threads 10  -data $TMP  -variance linExp  -optimize  -subgraph_iter_max 2  \
    -hybridness_min 1.2  -hybrid_parent_pairs $TMP.Saccharomyces2.hybrid_parent_pairs  -delete_hybrids $TMP.Saccharomyces2.hybrid  -dissim_boundary 0.675 \
    -delete_criterion_outliers $TMP.Saccharomyces2.criterion_outliers  -criterion_outlier_num_max 1 \
    -delete_deformation_outliers Saccharomyces2.deformation_outliers  -deformation_outlier_num_max 1 \
    -output_tree $TMP.tree > $TMP.out
  $THIS/../tsv/tsv_cut.sh $TMP.Saccharomyces2.hybrid "10 --complement" 0 | sort > $TMP.Saccharomyces2.hybrid.cut
  diff $TMP.Saccharomyces2.hybrid.cut $DATA/Saccharomyces2.hybrid
  $THIS/hybrid2list.sh $TMP.Saccharomyces2.hybrid > Saccharomyces2.hybrid.list
  $THIS/../tsv/tsv_cut.sh $TMP.Saccharomyces2.hybrid_parent_pairs "13 --complement" 0 | sort > $TMP.Saccharomyces2.hybrid_parent_pairs.cut
  diff $TMP.Saccharomyces2.hybrid_parent_pairs.cut $DATA/Saccharomyces2.hybrid_parent_pairs
  diff $TMP.Saccharomyces2.criterion_outliers $DATA/Saccharomyces2.criterion_outliers
  cut -f 1       Saccharomyces2.deformation_outliers > $TMP.Saccharomyces2.deformation_outliers
  cut -f 1 $DATA/Saccharomyces2.deformation_outliers > $TMP.data.Saccharomyces2.deformation_outliers
  diff $TMP.Saccharomyces2.deformation_outliers $TMP.data.Saccharomyces2.deformation_outliers
  rm Saccharomyces2.deformation_outliers 
  # Saccharomyces2.distTree
  $THIS/tree2obj.sh $TMP.tree > $TMP.list
  $THIS/../dm/dm2subset $TMP $TMP.list > $TMP.subset.dm
  $THIS/makeDistTree  -threads 10  -data $TMP.subset  -input_tree $TMP.tree  -variance linExp  -optimize  -reinsert  -subgraph_iter_max 10  > Saccharomyces2.distTree
  $THIS/distTree_compare_criteria.sh Saccharomyces2.distTree $DATA/Saccharomyces2.distTree
  rm Saccharomyces2.distTree

  set +o errexit
  N=`diff -y --suppress-common-lines Saccharomyces.hybrid.list Saccharomyces2.hybrid.list | wc -l`
  set -o errexit
  if [ $N -gt 1 ]; then  # ??
    diff Saccharomyces.hybrid.list Saccharomyces2.hybrid.list
  fi
  rm Saccharomyces.hybrid.list Saccharomyces2.hybrid.list
  diff $DATA/Saccharomyces2.criterion_outliers $DATA/Saccharomyces.criterion_outliers
  cut -f 1 $DATA/Saccharomyces2.deformation_outliers > $TMP.Saccharomyces2.deformation_outliers
  cut -f 1 $DATA/Saccharomyces.deformation_outliers  > $TMP.Saccharomyces.deformation_outliers
  diff $TMP.Saccharomyces2.deformation_outliers $TMP.Saccharomyces.deformation_outliers


  super_section "Many dissimilarity types"

  section "Virus9"
  gunzip -c $DATA/Virus9.dm.gz > $TMP.dm
  $THIS/makeDistTree  -qc  -data $TMP  -optimize  -subgraph_iter_max 100  -variance linExp  -variance_dissim  -output_dissim_coeff Virus9.coeff  -output_data Virus9-out  -output_tree Virus9.tree  1> $TMP.1 2> $TMP.out
  diff Virus9.coeff $DATA/Virus9.coeff
  rm Virus9.coeff
  A=`grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1`
  #
  $THIS/makeDistTree  -qc  -data Virus9-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Virus9-out.tree  1> $TMP.2 2> $TMP.out
  B=`grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1`
  if [ "$A" != "$B" ]; then
    error "$A != $B"
  fi
  $THIS/printDistTree  -qc  Virus9.tree      -order  -decimals 3  > Virus9.nw
  $THIS/printDistTree  -qc  Virus9-out.tree  -order  -decimals 3  > Virus9-out.nw
  diff Virus9.nw Virus9-out.nw
  rm Virus9.nw Virus9-out.nw
  rm Virus9.tree Virus9-out.tree
  rm Virus9-out.dm


  section "Virus110"
  gunzip -c $DATA/Virus110.dm.gz > $TMP.dm
  $THIS/makeDistTree  -data $TMP  -variance_min 0.005  -variance linExp  -variance_dissim  -optimize  -subgraph_iter_max 100  -output_data Virus110-out  -output_tree Virus110.tree 1> $TMP.1 
  A=`grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1`
  #
  $THIS/makeDistTree  -qc  -input_tree Virus110.tree  -data Virus110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Virus110-out1.tree 1> $TMP.2 2> $TMP.out
  B=`grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1`
  if [ "$A" != "$B" ]; then
    error "$A != $B"
  fi
  $THIS/printDistTree  -qc  Virus110.tree       -order  -decimals 1 | sed 's/,(/,\n(/g' > Virus110.nw
  $THIS/printDistTree  -qc  Virus110-out1.tree  -order  -decimals 1 | sed 's/,(/,\n(/g' > Virus110-out1.nw
  diff Virus110.nw Virus110-out1.nw || true
  #
  if false; then  # ??
    makeDistTree  -qc  -data Virus110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Virus110-out2.tree 1> $TMP.2 2> $TMP.out
    B=`grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1`
    if [ "$A" != "$B" ]; then
      error "$A != $B"
    fi
    $THIS/printDistTree  -qc  Virus110.tree       -order  -decimals 2  > Virus110.nw
    $THIS/printDistTree  -qc  Virus110-out2.tree  -order  -decimals 2  > Virus110-out2.nw
    diff Virus110.nw Virus110-out.nw
  fi
  #
  rm -f Virus110.nw Virus110-out1.nw Virus110-out2.nw
  rm -f Virus110.tree Virus110-out1.tree Virus110-out2.tree
  rm Virus110-out.dm
fi


rm -r $TMP*


success
