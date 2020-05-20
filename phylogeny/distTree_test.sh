#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: 69 min."
  exit 1
fi


TMP=`mktemp`
if [ ${AT_NCBI:-0} == 1 ]; then
  echo $TMP
fi


#if [ 1 == 0 ]; then  
echo ""
echo "mdsTree: Enterobacteriaceae ..."
rm -rf data/Enterobacteriaceae.dir/
$THIS/../dm/mdsTree.sh data/Enterobacteriaceae Conservation 2  &> /dev/null
$THIS/makeDistTree  -qc -input_tree data/Enterobacteriaceae.dir/  -data data/Enterobacteriaceae  -dissim_attr Conservation  -variance linExp  -optimize  -output_tree Enterobacteriaceae.tree > Enterobacteriaceae.distTree
$THIS/distTree_compare_criteria.sh Enterobacteriaceae.distTree data/Enterobacteriaceae.distTree
rm Enterobacteriaceae.distTree
rm -r data/Enterobacteriaceae.dir/

echo ""
echo "To Newick ..."
$THIS/printDistTree  -qc  -data data/Enterobacteriaceae  -dissim_attr Conservation  -variance linExp  Enterobacteriaceae.tree  -order  -decimals 4  -ext_name > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw data/Enterobacteriaceae.nw
rm Enterobacteriaceae.nw

rm Enterobacteriaceae.tree

echo ""
echo "From Newick ..."
$THIS/newick2tree -qc data/Enterobacteriaceae.nw > Enterobacteriaceae.tree
$THIS/printDistTree -qc Enterobacteriaceae.tree  -decimals 4  -ext_name > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw data/Enterobacteriaceae.nw
rm Enterobacteriaceae.nw
rm Enterobacteriaceae.tree

echo ""
echo "Perfect tree ..."
$THIS/makeDistTree  -qc  -data data/tree4  -variance linExp  -optimize  -output_tree tree4 | grep -v '^CHRON: ' > tree4.makeDistTree
$THIS/distTree_compare_criteria.sh tree4.makeDistTree data/tree4.makeDistTree
rm tree4.makeDistTree
$THIS/statDistTree tree4 > tree4.stat
rm tree4
diff tree4.stat data/tree4.stat
rm tree4.stat
echo "Verbose ..."
$THIS/makeDistTree -qc  -data data/tree4  -variance linExp  -optimize  -verbose 2 &> /dev/null

echo ""
echo "mdsTree: Random tree ..."
$THIS/makeDistTree  -qc   -data data/randomTree  -variance lin  -output_tree random-output.tree > /dev/null
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data data/randomTree  -variance lin  -optimize | grep -v '^CHRON: ' > randomTree.makeDistTree
$THIS/distTree_compare_criteria.sh randomTree.makeDistTree data/randomTree.makeDistTree
echo "Verbose ..."
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data data/randomTree  -variance lin  -verbose 2 &> /dev/null
rm randomTree.makeDistTree
rm random-output.tree

echo ""
echo "Salmonella ..."
# Check time ??
$THIS/makeDistTree  -qc  -data data/Salmonella  -variance linExp  -optimize  -subgraph_iter_max 10 \
  -delete_criterion_outliers Salmonella.criterion_outliers \
  -delete_deformation_outliers Salmonella.deformation_outliers \
  -delete_hybrids Salmonella.hybrids \
  | grep -v '^CHRON: ' > Salmonella.distTree
diff Salmonella.criterion_outliers data/Salmonella.criterion_outliers
rm Salmonella.criterion_outliers
$THIS/../sort.sh Salmonella.deformation_outliers
diff Salmonella.deformation_outliers data/Salmonella.deformation_outliers
rm Salmonella.deformation_outliers
diff Salmonella.hybrids data/Salmonella.hybrids
rm Salmonella.hybrids
$THIS/distTree_compare_criteria.sh Salmonella.distTree data/Salmonella.distTree
rm Salmonella.distTree

echo ""
echo "testDistTree ..."
$THIS/makeDistTree  -data data/Salmonella  -variance linExp  -optimize  -output_tree $TMP.tree > /dev/null
$THIS/testDistTree -qc  data/Salmonella  -input_tree  $TMP.tree  -variance linExp 

echo ""
echo "Salmonella: delete ..."
$THIS/makeDistTree  -qc  -data data/Salmonella  -variance linExp  -delete data/delete.list  -check_delete > /dev/null

if [ -d data/inc.ITS ]; then
  echo ""
  echo "ITS threads ..."
  $THIS/makeDistTree  -threads 10  -data data/inc.ITS/  -variance linExp  -variance_dissim  -optimize  -skip_len  -reinsert  -subgraph_iter_max 1  -noqual > ITS.distTree
  $THIS/distTree_compare_criteria.sh ITS.distTree data/ITS.distTree
  rm ITS.distTree
fi

echo ""
echo "Saccharomyces hybrids ..."
$THIS/makeDistTree -qc  -threads 3  -data data/Saccharomyces  -variance linExp  -optimize  -subgraph_iter_max 2  \
  -hybridness_min 1.2  -hybrid_parent_pairs Saccharomyces.hybrid_parent_pairs  -delete_hybrids Saccharomyces.hybrid  -dissim_boundary 0.675 \
  -delete_criterion_outliers Saccharomyces.criterion_outliers  -criterion_outlier_num_max 1 \
  -delete_deformation_outliers Saccharomyces.deformation_outliers  -deformation_outlier_num_max 1 \
  -output_tree $TMP.tree > /dev/null
diff Saccharomyces.hybrid data/Saccharomyces.hybrid
$THIS/hybrid2list.sh Saccharomyces.hybrid > Saccharomyces.hybrid.list
rm Saccharomyces.hybrid
diff Saccharomyces.hybrid_parent_pairs data/Saccharomyces.hybrid_parent_pairs
rm Saccharomyces.hybrid_parent_pairs
diff Saccharomyces.criterion_outliers data/Saccharomyces.criterion_outliers
rm Saccharomyces.criterion_outliers
diff Saccharomyces.deformation_outliers data/Saccharomyces.deformation_outliers
rm Saccharomyces.deformation_outliers
# Saccharomyces.distTree
$THIS/tree2obj.sh $TMP.tree > $TMP.list
$THIS/../dm/dm2subset data/Saccharomyces $TMP.list > $TMP.dm
$THIS/makeDistTree  -threads 3  -data $TMP  -input_tree $TMP.tree  -variance linExp  -optimize  -reinsert  -subgraph_iter_max 10  > Saccharomyces.distTree
$THIS/distTree_compare_criteria.sh Saccharomyces.distTree data/Saccharomyces.distTree
rm Saccharomyces.distTree

echo ""
echo "-variance_min ..."
# 0.0005 = average arc length / 100
$THIS/makeDistTree  -qc  -data data/Salmonella  -variance linExp  -variance_min 0.0005  -optimize  -subgraph_iter_max 10 \
  -delete_criterion_outliers Salmonella-var_min.criterion_outliers \
  -delete_deformation_outliers Salmonella-var_min.deformation_outliers \
  -delete_hybrids Salmonella-var_min.hybrids \
  | grep -v '^CHRON: ' > Salmonella-var_min.distTree
diff Salmonella-var_min.criterion_outliers data/Salmonella-var_min.criterion_outliers
rm Salmonella-var_min.criterion_outliers
diff Salmonella-var_min.deformation_outliers data/Salmonella-var_min.deformation_outliers
rm Salmonella-var_min.deformation_outliers
diff Salmonella-var_min.hybrids data/Salmonella-var_min.hybrids
rm Salmonella-var_min.hybrids
$THIS/distTree_compare_criteria.sh Salmonella-var_min.distTree data/Salmonella-var_min.distTree
rm Salmonella-var_min.distTree


if [ ${AT_NCBI:-0} == 0 ]; then
  rm Saccharomyces.hybrid.list
  rm $TMP*
  exit 0
fi


echo ""
echo ""
echo "Two dissimilarity types ..."

echo ""
echo "Salmonella2 ..."
$THIS/makeDistTree  -qc  -data data/Salmonella2  -variance linExp  -optimize  -subgraph_iter_max 10  \
  -output_dissim_coeff Salmonella2.dissim_coeff \
  -delete_criterion_outliers Salmonella2.criterion_outliers \
  -delete_deformation_outliers Salmonella2.deformation_outliers \
  -delete_hybrids Salmonella2.hybrids \
  | grep -v '^CHRON: ' > Salmonella2.distTree
diff Salmonella2.dissim_coeff data/Salmonella2.dissim_coeff
rm Salmonella2.dissim_coeff
diff Salmonella2.criterion_outliers data/Salmonella2.criterion_outliers
rm Salmonella2.criterion_outliers
$THIS/../sort.sh Salmonella2.deformation_outliers
diff Salmonella2.deformation_outliers data/Salmonella2.deformation_outliers
rm Salmonella2.deformation_outliers
diff Salmonella2.hybrids data/Salmonella2.hybrids
rm Salmonella2.hybrids
$THIS/distTree_compare_criteria.sh Salmonella2.distTree data/Salmonella2.distTree
rm Salmonella2.distTree

echo ""
echo "Saccharomyces hybrids ..."
$THIS/makeDistTree -qc  -threads 3  -data data/Saccharomyces2  -variance linExp  -optimize  -subgraph_iter_max 2  \
  -hybridness_min 1.2  -hybrid_parent_pairs Saccharomyces2.hybrid_parent_pairs  -delete_hybrids Saccharomyces2.hybrid  -dissim_boundary 0.675 \
  -delete_criterion_outliers Saccharomyces2.criterion_outliers  -criterion_outlier_num_max 1 \
  -delete_deformation_outliers Saccharomyces2.deformation_outliers  -deformation_outlier_num_max 1 \
  -output_tree $TMP.tree > /dev/null
$THIS/../rm_col.sh Saccharomyces2.hybrid 10
$THIS/../sort.sh Saccharomyces2.hybrid
diff Saccharomyces2.hybrid data/Saccharomyces2.hybrid
$THIS/hybrid2list.sh Saccharomyces2.hybrid > Saccharomyces2.hybrid.list
rm Saccharomyces2.hybrid
$THIS/../rm_col.sh Saccharomyces2.hybrid_parent_pairs 13
$THIS/../sort.sh Saccharomyces2.hybrid_parent_pairs
diff Saccharomyces2.hybrid_parent_pairs data/Saccharomyces2.hybrid_parent_pairs
rm Saccharomyces2.hybrid_parent_pairs
diff Saccharomyces2.criterion_outliers data/Saccharomyces2.criterion_outliers
rm Saccharomyces2.criterion_outliers
diff Saccharomyces2.deformation_outliers data/Saccharomyces2.deformation_outliers
rm Saccharomyces2.deformation_outliers
# Saccharomyces2.distTree
$THIS/tree2obj.sh $TMP.tree > $TMP.list
$THIS/../dm/dm2subset data/Saccharomyces2 $TMP.list > $TMP.dm
$THIS/makeDistTree  -threads 3  -data $TMP  -input_tree $TMP.tree  -variance linExp  -optimize  -reinsert  -subgraph_iter_max 10  > Saccharomyces2.distTree
$THIS/distTree_compare_criteria.sh Saccharomyces2.distTree data/Saccharomyces2.distTree
rm Saccharomyces2.distTree

set +o errexit
N=`diff -y --suppress-common-lines Saccharomyces.hybrid.list Saccharomyces2.hybrid.list | wc -l`
set -o errexit
if [ $N -gt 1 ]; then  # ??
  diff Saccharomyces.hybrid.list Saccharomyces2.hybrid.list
fi
rm Saccharomyces.hybrid.list Saccharomyces2.hybrid.list
diff data/Saccharomyces2.criterion_outliers data/Saccharomyces.criterion_outliers
diff data/Saccharomyces2.deformation_outliers data/Saccharomyces.deformation_outliers


echo ""
echo ""
echo "Many dissimilarity types ..."

echo ""
echo "Virus9 ..."
$THIS/makeDistTree  -qc  -data data/Virus9  -optimize  -variance linExp  -variance_dissim  -output_dissim_coeff Virus9.coeff  -output_data Virus9-out  -output_tree Virus9.tree  1> $TMP.1 2> /dev/null
diff Virus9.coeff data/Virus9.coeff
rm Virus9.coeff
A=`grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1`
#
$THIS/makeDistTree  -qc  -data Virus9-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Virus9-out.tree  1> $TMP.2 2> /dev/null
B=`grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1`
if [ "$A" != "$B" ]; then
  echo "$A"
  echo "$B"
  exit 1
fi
$THIS/printDistTree  -qc  Virus9.tree      -order  -decimals 3  > Virus9.nw
$THIS/printDistTree  -qc  Virus9-out.tree  -order  -decimals 3  > Virus9-out.nw
diff Virus9.nw Virus9-out.nw
rm Virus9.nw Virus9-out.nw
rm Virus9.tree Virus9-out.tree
rm Virus9-out.dm


echo ""
echo "Virus110 ..."
$THIS/makeDistTree  -data data/Virus110  -variance_min 0.005  -variance linExp  -variance_dissim  -optimize  -subgraph_iter_max 100  -output_data Virus110-out  -output_tree Virus110.tree 1> $TMP.1 
A=`grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1`
#
$THIS/makeDistTree  -qc  -input_tree Virus110.tree  -data Virus110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Virus110-out1.tree 1> $TMP.2 2> /dev/null
B=`grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1`
if [ "$A" != "$B" ]; then
  echo "$A"
  echo "$B"
  exit 1
fi
$THIS/printDistTree  -qc  Virus110.tree       -order  -decimals 1 | sed -e 's/,(/,\n(/g' > Virus110.nw
$THIS/printDistTree  -qc  Virus110-out1.tree  -order  -decimals 1 | sed -e 's/,(/,\n(/g' > Virus110-out1.nw
diff Virus110.nw Virus110-out1.nw || true
#
if [ 0 == 1 ]; then  # ??
  makeDistTree  -qc  -data Virus110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Virus110-out2.tree 1> $TMP.2 2> /dev/null
  B=`grep -w '^OUTPUT' -A 1 $TMP.2 | tail -1`
  if [ "$A" != "$B" ]; then
    echo "$A"
    echo "$B"
    exit 1
  fi
  $THIS/printDistTree  -qc  Virus110.tree       -order  -decimals 2  > Virus110.nw
  $THIS/printDistTree  -qc  Virus110-out2.tree  -order  -decimals 2  > Virus110-out2.nw
  diff Virus110.nw Virus110-out.nw
fi
#
rm -f Virus110.nw Virus110-out1.nw Virus110-out2.nw
rm -f Virus110.tree Virus110-out1.tree Virus110-out2.tree
rm Virus110-out.dm


rm $TMP*
