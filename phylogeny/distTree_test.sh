#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: ~60 min."
  exit 1
fi


TMP=`mktemp`
echo $TMP


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
$THIS/printDistTree  -qc  -data data/Enterobacteriaceae  -dissim_attr Conservation  -variance linExp  Enterobacteriaceae.tree  -order  -decimals 4 > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw data/Enterobacteriaceae.nw
rm Enterobacteriaceae.nw

rm Enterobacteriaceae.tree

echo ""
echo "From Newick ..."
$THIS/newick2tree -qc data/Enterobacteriaceae.nw > Enterobacteriaceae.tree
$THIS/printDistTree -qc Enterobacteriaceae.tree  -decimals 4 > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw data/Enterobacteriaceae.nw
rm Enterobacteriaceae.nw
rm Enterobacteriaceae.tree

echo ""
echo "mdsTree: Mycobacterium_tuberculosis ..."
rm -rf data/Mycobacterium_tuberculosis.dir/
$THIS/../dm/mdsTree.sh data/Mycobacterium_tuberculosis ANI 2 &> /dev/null
rm -r data/Mycobacterium_tuberculosis.dir/

echo ""
echo "Perfect tree ..."
$THIS/makeDistTree  -qc  -data data/tree4  -variance linExp  -optimize  -output_tree tree4 | grep -v '^CHRON: ' > tree4.makeDistTree
diff tree4.makeDistTree data/tree4.makeDistTree
rm tree4.makeDistTree
$THIS/statDistTree tree4 > tree4.stat
rm tree4
diff tree4.stat data/tree4.stat
rm tree4.stat
echo "Verbose ..."
$THIS/makeDistTree -qc  -data data/tree4  -variance linExp  -optimize  -verbose 2 &> /dev/null

echo ""
echo "mdsTree: Random tree ..."
rm -rf data/randomTree.dir/
$THIS/../dm/mdsTree.sh data/randomTree dist 2 &> /dev/null
$THIS/makeDistTree  -qc  -input_tree data/randomTree.dir/  -data data/randomTree  -variance lin  -output_tree random-output.tree > /dev/null
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data data/randomTree  -variance lin  -optimize | grep -v '^CHRON: ' > randomTree.makeDistTree
diff randomTree.makeDistTree data/randomTree.makeDistTree
echo "Verbose ..."
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data data/randomTree  -variance lin  -verbose 2 &> /dev/null
rm -r data/randomTree.dir/
rm randomTree.makeDistTree
rm random-output.tree

echo ""
echo "prot-identical_comm: subgraphs ..."
# Check time ??
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -variance linExp  -optimize  \
  -delete_outliers prot-identical_comm.outliers \
  -delete_hybrids prot-identical_comm.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm.distTree
diff prot-identical_comm.outliers data/prot-identical_comm.outliers
rm prot-identical_comm.outliers
diff prot-identical_comm.hybrids data/prot-identical_comm.hybrids
rm prot-identical_comm.hybrids
$THIS/distTree_compare_criteria.sh prot-identical_comm.distTree data/prot-identical_comm.distTree
rm prot-identical_comm.distTree

echo ""
echo "prot-identical_comm: subgraphs, delete ..."
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -variance linExp  -delete data/delete.list  -check_delete > /dev/null

echo ""
echo "ITS threads ..."
# -qc: 40 min.
$THIS/makeDistTree  -data data/inc.ITS/  -variance linExp  -optimize  -skip_len  -subgraph_iter_max 1  -noqual  -threads 5 > ITS.distTree
$THIS/distTree_compare_criteria.sh ITS.distTree data/ITS.distTree
rm ITS.distTree


echo ""
echo "Saccharomyces hybrids ..."
$THIS/makeDistTree -qc  -threads 3  -data data/Saccharomyces  -variance linExp  -optimize  -subgraph_iter_max 2  \
  -hybridness_min 1.2  -hybrid_parent_pairs Saccharomyces.hybrid_parent_pairs  -delete_hybrids Saccharomyces.hybrid  -dissim_boundary 0.675 \
  -delete_outliers Saccharomyces.outliers  -outlier_num_max 1 \
  > Saccharomyces.distTree
diff Saccharomyces.hybrid data/Saccharomyces.hybrid
$THIS/hybrid2list.sh Saccharomyces.hybrid > Saccharomyces.hybrid.list
rm Saccharomyces.hybrid
diff Saccharomyces.hybrid_parent_pairs data/Saccharomyces.hybrid_parent_pairs
rm Saccharomyces.hybrid_parent_pairs
diff Saccharomyces.outliers data/Saccharomyces.outliers
rm Saccharomyces.outliers
$THIS/distTree_compare_criteria.sh Saccharomyces.distTree data/Saccharomyces.distTree
rm Saccharomyces.distTree


echo ""
echo "-var_min ..."
# 0.0005 = average arc length / 100
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -variance linExp  -var_min 0.0005  -optimize  \
  -delete_outliers prot-identical_comm-var_min.outliers \
  -delete_hybrids prot-identical_comm-var_min.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm-var_min.distTree
diff prot-identical_comm-var_min.outliers data/prot-identical_comm-var_min.outliers
rm prot-identical_comm-var_min.outliers
diff prot-identical_comm-var_min.hybrids data/prot-identical_comm-var_min.hybrids
rm prot-identical_comm-var_min.hybrids
$THIS/distTree_compare_criteria.sh prot-identical_comm-var_min.distTree data/prot-identical_comm-var_min.distTree
rm prot-identical_comm-var_min.distTree


echo ""
echo ""
echo "Two dissimilarity types ..."

echo ""
echo "prot-identical_comm: subgraphs ..."
makeDistTree  -qc  -data data/prot-identical_comm2  -variance linExp  -optimize  -subgraph_iter_max 10  \
  -output_dissim_coeff prot-identical_comm2.dissim_coeff \
  -delete_outliers prot-identical_comm2.outliers \
  -delete_hybrids prot-identical_comm2.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm2.distTree
diff prot-identical_comm2.dissim_coeff data/prot-identical_comm2.dissim_coeff
rm prot-identical_comm2.dissim_coeff
diff prot-identical_comm2.outliers data/prot-identical_comm2.outliers
rm prot-identical_comm2.outliers
diff prot-identical_comm2.hybrids data/prot-identical_comm2.hybrids
rm prot-identical_comm2.hybrids
$THIS/distTree_compare_criteria.sh prot-identical_comm2.distTree data/prot-identical_comm2.distTree
rm prot-identical_comm2.distTree

echo ""
echo "Saccharomyces hybrids ..."
$THIS/makeDistTree -qc  -threads 3  -data data/Saccharomyces2  -variance linExp  -optimize  -subgraph_iter_max 2  \
  -hybridness_min 1.2  -hybrid_parent_pairs Saccharomyces2.hybrid_parent_pairs  -delete_hybrids Saccharomyces2.hybrid  -dissim_boundary 0.675 \
  -delete_outliers Saccharomyces2.outliers  -outlier_num_max 1 \
  > Saccharomyces2.distTree
diff Saccharomyces2.hybrid data/Saccharomyces2.hybrid
$THIS/hybrid2list.sh Saccharomyces2.hybrid > Saccharomyces2.hybrid.list
rm Saccharomyces2.hybrid
diff Saccharomyces2.hybrid_parent_pairs data/Saccharomyces2.hybrid_parent_pairs
rm Saccharomyces2.hybrid_parent_pairs
diff Saccharomyces2.outliers data/Saccharomyces2.outliers
rm Saccharomyces2.outliers
$THIS/distTree_compare_criteria.sh Saccharomyces2.distTree data/Saccharomyces2.distTree
rm Saccharomyces2.distTree

N=`diff -y --suppress-common-lines Saccharomyces.hybrid.list Saccharomyces2.hybrid.list | wc -l`
if [ $N -gt 1 ]; then
  diff Saccharomyces.hybrid.list Saccharomyces2.hybrid.list
fi
rm Saccharomyces.hybrid.list Saccharomyces2.hybrid.list
diff data/Saccharomyces2.outliers data/Saccharomyces.outliers


echo ""
echo ""
echo "Many dissimilarity types ..."

echo ""
echo "Wolf9 ..."
makeDistTree  -qc  -data data/Wolf9  -optimize  -variance linExp  -output_dissim_coeff Wolf9.coeff  -output_data Wolf9-out  -output_tree Wolf9.tree 1> $TMP.1 2> /dev/null
diff Wolf9.coeff data/Wolf9.coeff
rm Wolf9.coeff
A=`grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1`
makeDistTree  -qc  -data Wolf9-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Wolf9-out.tree 1> $TMP.2 2> /dev/null
B=`grep -w ^OUTPUT -A 1 $TMP.2 | tail -1`
if [ "$A" != "$B" ]; then
  echo "$A"
  echo "$B"
  exit 1
fi
printDistTree  -qc  Wolf9.tree      -order  -decimals 3  -min_name > Wolf9.nw
printDistTree  -qc  Wolf9-out.tree  -order  -decimals 3  -min_name > Wolf9-out.nw
diff Wolf9.nw Wolf9-out.nw
rm Wolf9.nw Wolf9-out.nw
rm Wolf9.tree Wolf9-out.tree
rm Wolf9-out.dm

echo ""
echo "Wolf110 ..."
# 643 iterations, 14 min.
makeDistTree  -qc  -data data/Wolf110  -var_min 0.005  -variance linExp  -optimize  -output_data Wolf110-out  -output_tree Wolf110.tree 1> $TMP.1 
A=`grep -w '^Error between dissimilarities' -A 1 $TMP.1 | tail -1`
#
makeDistTree  -qc  -input_tree Wolf110.tree  -data Wolf110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Wolf110-out.tree 1> $TMP.2 2> /dev/null
B=`grep -w ^OUTPUT -A 1 $TMP.2 | tail -1`
if [ "$A" != "$B" ]; then
  echo "$A"
  echo "$B"
  exit 1
fi
printDistTree  -qc  Wolf110.tree      -order  -decimals 2  -min_name > Wolf110.nw
printDistTree  -qc  Wolf110-out.tree  -order  -decimals 2  -min_name > Wolf110-out.nw
diff Wolf110.nw Wolf110-out.nw
rm Wolf110.nw Wolf110-out.nw
rm Wolf110.tree Wolf110-out.tree
rm Wolf110-out.dm
#
# makeDistTree  -qc  -data Wolf110-out  -dissim_attr dissim  -weight_attr weight  -optimize  -output_tree Wolf110-out.tree 1> $TMP.2 2> /dev/null
# ... ??



if [ 1 == 0 ]; then
	echo ""
	echo "prot-identical_comm: whole ..."
	# Time: 7 min.
	$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -optimize  -whole  \
	  -delete_outliers prot-identical_comm.outliers-whole \
	  | grep -v '^CHRON: ' > prot-identical_comm.distTree-whole
	diff prot-identical_comm.outliers-whole data/prot-identical_comm.outliers-whole
	rm prot-identical_comm.outliers-whole
	diff prot-identical_comm.distTree-whole data/prot-identical_comm.distTree-whole
	rm prot-identical_comm.distTree-whole
fi


rm $TMP*
