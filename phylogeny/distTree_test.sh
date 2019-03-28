#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: 31 min."
  exit 1
fi


TMP=`mktemp`
echo $TMP


#if [ 1 == 0 ]; then  
echo ""
echo "mdsTree: Enterobacteriaceae ..."
rm -rf data/Enterobacteriaceae.dir/
$THIS/../dm/mdsTree.sh data/Enterobacteriaceae Conservation 2 &> /dev/null
$THIS/makeDistTree  -qc -input_tree data/Enterobacteriaceae.dir/  -data data/Enterobacteriaceae  -dissim_attr Conservation  -optimize  -output_tree Enterobacteriaceae.tree > /dev/null
rm -r data/Enterobacteriaceae.dir/

echo ""
echo "To Newick ..."
$THIS/printDistTree  -qc  -data data/Enterobacteriaceae  -dissim_attr Conservation  Enterobacteriaceae.tree  -order  -decimals 4 > Enterobacteriaceae.nw
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
$THIS/makeDistTree  -qc  -data data/tree4  -optimize | grep -v '^CHRON: ' > tree4.makeDistTree
diff tree4.makeDistTree data/tree4.makeDistTree
rm tree4.makeDistTree
echo "Verbose ..."
$THIS/makeDistTree -qc  -data data/tree4  -optimize  -verbose 2 &> /dev/null

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
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -optimize  \
  -delete_outliers prot-identical_comm.outliers \
  -delete_hybrids prot-identical_comm.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm.distTree
diff prot-identical_comm.outliers data/prot-identical_comm.outliers
rm prot-identical_comm.outliers
diff prot-identical_comm.hybrids data/prot-identical_comm.hybrids
rm prot-identical_comm.hybrids
grep OUTPUT -A 1      prot-identical_comm.distTree > $TMP.1
grep OUTPUT -A 1 data/prot-identical_comm.distTree > $TMP.2
diff $TMP.1 $TMP.2
rm prot-identical_comm.distTree

echo ""
echo "prot-identical_comm: subgraphs, -variance sqr ..."
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -optimize  -threads 3  -variance sqr | grep -v '^CHRON: ' > prot-identical_comm-sqr.distTree
diff prot-identical_comm-sqr.distTree data/prot-identical_comm-sqr.distTree
rm prot-identical_comm-sqr.distTree

echo ""
echo "prot-identical_comm: subgraphs, delete ..."
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -delete data/delete.list  -check_delete > /dev/null

echo ""
echo "ITS threads ..."
$THIS/makeDistTree -data data/inc.ITS/  -variance lin  -optimize  -skip_len  -subgraph_iter_max 1  -noqual  -threads 5 > /dev/null

echo ""
echo "Saccharomyces hybrids ..."
$THIS/makeDistTree -qc  -threads 3  -data data/Saccharomyces  -optimize  -subgraph_iter_max 2  \
  -hybrid_parent_pairs Saccharomyces.hybrid_parent_pairs  -delete_hybrids Saccharomyces.hybrid  -dissim_boundary 0.675 \
  -delete_outliers Saccharomyces.outliers  -max_outlier_num 1 \
  > /dev/null
diff Saccharomyces.hybrid data/Saccharomyces.hybrid
rm Saccharomyces.hybrid
diff Saccharomyces.hybrid_parent_pairs data/Saccharomyces.hybrid_parent_pairs
cut -f 1,2,3,7,8,9,10,11,12 Saccharomyces.hybrid_parent_pairs > Saccharomyces.hybrid_parent_pairs.stable
rm Saccharomyces.hybrid_parent_pairs
diff Saccharomyces.outliers data/Saccharomyces.outliers
rm Saccharomyces.outliers


echo ""
echo "-min_var ..."
# 0.0005 = average arc length / 100
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -min_var 0.0005  -optimize  \
  -delete_outliers prot-identical_comm-min_var.outliers \
  -delete_hybrids prot-identical_comm-min_var.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm-min_var.distTree
diff prot-identical_comm-min_var.outliers data/prot-identical_comm-min_var.outliers
rm prot-identical_comm-min_var.outliers
diff prot-identical_comm-min_var.hybrids data/prot-identical_comm-min_var.hybrids
rm prot-identical_comm-min_var.hybrids
grep OUTPUT -A 1      prot-identical_comm-min_var.distTree > $TMP.1
grep OUTPUT -A 1 data/prot-identical_comm-min_var.distTree > $TMP.2
diff $TMP.1 $TMP.2
rm prot-identical_comm-min_var.distTree


echo ""
echo ""
echo "Two dissimilarity types ..."

echo ""
echo "prot-identical_comm: subgraphs ..."
makeDistTree  -qc  -data data/prot-identical_comm2  -optimize  -output_dissim_coeff prot-identical_comm2.dissim_coeff \
  -delete_outliers prot-identical_comm2.outliers \
  -delete_hybrids prot-identical_comm2.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm2.distTree
diff prot-identical_comm2.dissim_coeff data/prot-identical_comm2.dissim_coeff
rm prot-identical_comm2.dissim_coeff
diff prot-identical_comm2.outliers data/prot-identical_comm2.outliers
rm prot-identical_comm2.outliers
diff prot-identical_comm2.hybrids data/prot-identical_comm2.hybrids
rm prot-identical_comm2.hybrids
grep OUTPUT -A 1      prot-identical_comm2.distTree > $TMP.1
grep OUTPUT -A 1 data/prot-identical_comm2.distTree > $TMP.2
diff $TMP.1 $TMP.2
rm prot-identical_comm2.distTree

echo ""
echo "Saccharomyces hybrids ..."
$THIS/makeDistTree -qc  -threads 3  -data data/Saccharomyces2  -optimize  -subgraph_iter_max 2  \
  -hybrid_parent_pairs Saccharomyces2.hybrid_parent_pairs  -delete_hybrids Saccharomyces2.hybrid  -dissim_boundary 0.675 \
  -delete_outliers Saccharomyces2.outliers  -max_outlier_num 1 \
  > /dev/null
diff Saccharomyces2.hybrid data/Saccharomyces2.hybrid
rm Saccharomyces2.hybrid
diff Saccharomyces2.hybrid_parent_pairs data/Saccharomyces2.hybrid_parent_pairs
cut -f 1,2,3,7,8,9,10,11,12 Saccharomyces2.hybrid_parent_pairs > Saccharomyces2.hybrid_parent_pairs.stable
rm Saccharomyces2.hybrid_parent_pairs
diff Saccharomyces2.outliers data/Saccharomyces2.outliers
rm Saccharomyces2.outliers

diff Saccharomyces.hybrid_parent_pairs.stable Saccharomyces2.hybrid_parent_pairs.stable
rm Saccharomyces.hybrid_parent_pairs.stable Saccharomyces2.hybrid_parent_pairs.stable
diff data/Saccharomyces2.outliers data/Saccharomyces.outliers


echo ""
echo ""
echo "Many dissimilarity types ..."
makeDistTree  -qc  -data data/Wolf9  -optimize  -output_dissim_coeff Wolf9.coeff  -output_data Wolf9-out  -output_tree Wolf9.tree 1> $TMP.1 2> /dev/null
diff Wolf9.coeff data/Wolf9.coeff
rm Wolf9.coeff
A=`grep -w ^Unoptimizable -A 1 $TMP.1 | tail -1`
makeDistTree  -qc  -data Wolf9-out  -dissim_attr dissim  -var_attr var  -optimize  -output_tree Wolf9-out.tree 1> $TMP.2 2> /dev/null
B=`grep -w ^OUTPUT -A 1 $TMP.2 | tail -1`
if [ "$A" != "$B" ]; then
  echo "$A"
  echo "$B"
  exit 1
fi
printDistTree  -qc  Wolf9.tree      -order  -decimals 4  -min_name > Wolf9.nw
printDistTree  -qc  Wolf9-out.tree  -order  -decimals 4  -min_name > Wolf9-out.nw
diff Wolf9.nw Wolf9-out.nw
rm Wolf9.nw Wolf9-out.nw
rm Wolf9.tree Wolf9-out.tree
rm Wolf9-out.dm



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
