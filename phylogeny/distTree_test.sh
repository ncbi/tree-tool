#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  exit 1
fi


#if [ 1 == 0 ]; then  
echo ""
echo "mdsTree: Enterobacteriaceae ..."
rm -rf data/Enterobacteriaceae.dir/
mdsTree.sh data/Enterobacteriaceae Conservation 2 &> /dev/null

makeDistTree  -qc -input_tree data/Enterobacteriaceae.dir/  -data data/Enterobacteriaceae  -dissim Conservation  -optimize  -output_tree Enterobacteriaceae.tree > /dev/null
rm -r data/Enterobacteriaceae.dir/

printDistTree  -qc  -data data/Enterobacteriaceae  -dissim Conservation Enterobacteriaceae.tree > Enterobacteriaceae.nw
diff Enterobacteriaceae.nw data/Enterobacteriaceae.nw
rm Enterobacteriaceae.tree
rm Enterobacteriaceae.nw


echo ""
echo "mdsTree: Mycobacterium_tuberculosis ..."
rm -r -f data/Mycobacterium_tuberculosis.dir/
mdsTree.sh data/Mycobacterium_tuberculosis ANI 2 &> /dev/null
rm -r data/Mycobacterium_tuberculosis.dir/


echo ""
echo "Perfect tree ..."
makeDistTree -qc -data data/tree4 -dissim dist | grep -v '^CHRON: ' > tree4.makeDistTree
diff tree4.makeDistTree data/tree4.makeDistTree
rm tree4.makeDistTree
echo "Verbose..."
makeDistTree -qc -data data/tree4 -dissim dist  &> /dev/null


echo ""
echo "mdsTree: Random tree ..."
rm -r -f data/randomTree.dir/
mdsTree.sh data/randomTree dist 2 &> /dev/null
makeDistTree  -qc  -input_tree data/randomTree.dir/  -data data/randomTree  -dissim dist  -variance lin  -optimize  -output_tree random-output.tree > /dev/null
makeDistTree  -qc  -input_tree random-output.tree  -data data/randomTree  -dissim dist  -variance lin | grep -v '^CHRON: ' > randomTree.makeDistTree
diff randomTree.makeDistTree data/randomTree.makeDistTree
echo "Verbose..."
makeDistTree  -qc  -verbose 2  -input_tree random-output.tree  -data data/randomTree  -dissim dist  -variance lin  &> /dev/null
rm -r data/randomTree.dir/
rm randomTree.makeDistTree
rm random-output.tree


echo ""
echo "prot-identical_comm: subgraphs ..."
# Check time ??
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  \
  -optimize  \
  -delete_outliers prot-identical_comm.outliers \
  -delete_hybrids prot-identical_comm.hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm.distTree
diff prot-identical_comm.outliers data/prot-identical_comm.outliers
rm prot-identical_comm.outliers
diff prot-identical_comm.hybrids data/prot-identical_comm.hybrids
rm prot-identical_comm.hybrids
diff prot-identical_comm.distTree data/prot-identical_comm.distTree
rm prot-identical_comm.distTree

echo ""
echo "prot-identical_comm: subgraphs, threads ..."
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  -optimize  -threads 3 > /dev/null

echo ""
echo "prot-identical_comm: subgraphs, delete ..."
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  -delete data/delete.list > /dev/null

echo ""
echo "Saccharomyces hybrids ..."
makeDistTree -qc  -threads 3  -data data/Saccharomyces  -dissim cons  \
  -optimize  -subgraph_iter_max 2  -reinsert \
  -hybrid_parent_pairs Saccharomyces.hybrid_parent_pairs  -delete_hybrids Saccharomyces.hybrid  -delete_all_hybrids  -dissim_boundary 0.675 > /dev/null
diff Saccharomyces.hybrid data/Saccharomyces.hybrid
rm Saccharomyces.hybrid
diff Saccharomyces.hybrid_parent_pairs data/Saccharomyces.hybrid_parent_pairs
rm Saccharomyces.hybrid_parent_pairs

if [ 1 == 0 ]; then
	echo ""
	echo "prot-identical_comm: whole ..."
	# Time: 7 min.
	makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  \
	  -optimize  -whole  \
	  -delete_outliers prot-identical_comm.outliers-whole \
	  | grep -v '^CHRON: ' > prot-identical_comm.distTree-whole
	diff prot-identical_comm.outliers-whole data/prot-identical_comm.outliers-whole
	rm prot-identical_comm.outliers-whole
	diff prot-identical_comm.distTree-whole data/prot-identical_comm.distTree-whole
	rm prot-identical_comm.distTree-whole
fi

