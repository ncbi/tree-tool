#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: 14 min."
  exit 1
fi


#if [ 1 == 0 ]; then  
echo ""
echo "mdsTree: Enterobacteriaceae ..."
rm -rf data/Enterobacteriaceae.dir/
$THIS/../dm/mdsTree.sh data/Enterobacteriaceae Conservation 2 &> /dev/null
$THIS/makeDistTree  -qc -input_tree data/Enterobacteriaceae.dir/  -data data/Enterobacteriaceae  -dissim Conservation  -optimize  -output_tree Enterobacteriaceae.tree > /dev/null
rm -r data/Enterobacteriaceae.dir/

echo ""
echo "To Newick ..."
$THIS/printDistTree  -qc  -data data/Enterobacteriaceae  -dissim Conservation  Enterobacteriaceae.tree  -order > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw data/Enterobacteriaceae.nw
rm Enterobacteriaceae.tree
rm Enterobacteriaceae.nw

echo ""
echo "From Newick ..."
$THIS/newick2tree -qc data/Enterobacteriaceae.nw > Enterobacteriaceae.tree
$THIS/printDistTree -qc Enterobacteriaceae.tree > Enterobacteriaceae.nw
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
$THIS/makeDistTree  -qc  -data data/tree4  -dissim dist -optimize | grep -v '^CHRON: ' > tree4.makeDistTree
diff tree4.makeDistTree data/tree4.makeDistTree
rm tree4.makeDistTree
echo "Verbose ..."
$THIS/makeDistTree -qc  -data data/tree4  -dissim dist  -optimize  -verbose 2 &> /dev/null

echo ""
echo "mdsTree: Random tree ..."
rm -rf data/randomTree.dir/
$THIS/../dm/mdsTree.sh data/randomTree dist 2 &> /dev/null
$THIS/makeDistTree  -qc  -input_tree data/randomTree.dir/  -data data/randomTree  -dissim dist  -variance lin  -output_tree random-output.tree > /dev/null
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data data/randomTree  -dissim dist  -variance lin  -optimize | grep -v '^CHRON: ' > randomTree.makeDistTree
diff randomTree.makeDistTree data/randomTree.makeDistTree
echo "Verbose ..."
$THIS/makeDistTree  -qc  -input_tree random-output.tree    -data data/randomTree  -dissim dist  -variance lin  -verbose 2 &> /dev/null
rm -r data/randomTree.dir/
rm randomTree.makeDistTree
rm random-output.tree

echo ""
echo "prot-identical_comm: subgraphs ..."
# Check time ??
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  \
  -optimize  \
  -delete_outliers prot-identical_comm.outliers \
  -delete_hybrids prot-identical_comm.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm.distTree
diff prot-identical_comm.outliers data/prot-identical_comm.outliers
rm prot-identical_comm.outliers
diff prot-identical_comm.hybrids data/prot-identical_comm.hybrids
rm prot-identical_comm.hybrids
TMP=`mktemp`
grep OUTPUT -A 1      prot-identical_comm.distTree > $TMP.1
grep OUTPUT -A 1 data/prot-identical_comm.distTree > $TMP.2
diff $TMP.1 $TMP.2
rm $TMP*
rm prot-identical_comm.distTree

echo ""
echo "prot-identical_comm: subgraphs, -variance sqr ..."
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  -optimize  -threads 3  -variance sqr | grep -v '^CHRON: ' > prot-identical_comm-sqr.distTree
diff prot-identical_comm-sqr.distTree data/prot-identical_comm-sqr.distTree
rm prot-identical_comm-sqr.distTree

echo ""
echo "prot-identical_comm: subgraphs, delete ..."
$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  -delete data/delete.list  -check_delete > /dev/null

echo ""
echo "ITS threads ..."
$THIS/makeDistTree -data data/inc.ITS/  -variance lin  -optimize  -skip_len  -subgraph_iter_max 1  -noqual  -threads 5 > /dev/null

echo ""
echo "Saccharomyces hybrids ..."
$THIS/makeDistTree -qc  -threads 3  -data data/Saccharomyces  -dissim cons  \
  -optimize  -subgraph_iter_max 2  -reinsert \
  -hybrid_parent_pairs Saccharomyces.hybrid_parent_pairs  -delete_hybrids Saccharomyces.hybrid  -delete_all_hybrids  -dissim_boundary 0.675 \
  -delete_outliers Saccharomyces.outliers  -max_outlier_num 1 \
  > /dev/null
diff Saccharomyces.hybrid data/Saccharomyces.hybrid
rm Saccharomyces.hybrid
diff Saccharomyces.hybrid_parent_pairs data/Saccharomyces.hybrid_parent_pairs
rm Saccharomyces.hybrid_parent_pairs
diff Saccharomyces.outliers data/Saccharomyces.outliers
rm Saccharomyces.outliers


if [ 1 == 0 ]; then
	echo ""
	echo "prot-identical_comm: whole ..."
	# Time: 7 min.
	$THIS/makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  \
	  -optimize  -whole  \
	  -delete_outliers prot-identical_comm.outliers-whole \
	  | grep -v '^CHRON: ' > prot-identical_comm.distTree-whole
	diff prot-identical_comm.outliers-whole data/prot-identical_comm.outliers-whole
	rm prot-identical_comm.outliers-whole
	diff prot-identical_comm.distTree-whole data/prot-identical_comm.distTree-whole
	rm prot-identical_comm.distTree-whole
fi

