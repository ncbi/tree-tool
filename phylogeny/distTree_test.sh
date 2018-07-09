#!/bin/csh -f

if ($# != 1) then
  echo "#1: go"
  exit 1
endif


echo ""
echo "mdsTree: Enterobacteriaceae ..."
rm -r -f data/Enterobacteriaceae.dir/
mdsTree.sh data/Enterobacteriaceae Conservation 2 >& /dev/null
if ($?) exit 1

makeDistTree  -qc -input_tree data/Enterobacteriaceae.dir/  -data data/Enterobacteriaceae  -dissim Conservation  -optimize  -output_tree Enterobacteriaceae.tree > /dev/null
if ($?) exit 1
rm -r data/Enterobacteriaceae.dir/

printDistTree  -qc  -data data/Enterobacteriaceae  -dissim Conservation Enterobacteriaceae.tree > Enterobacteriaceae.nw
if ($?) exit 1
diff Enterobacteriaceae.nw data/Enterobacteriaceae.nw
if ($?) exit 1
rm Enterobacteriaceae.tree
rm Enterobacteriaceae.nw


echo ""
echo "mdsTree: Mycobacterium_tuberculosis ..."
rm -r -f data/Mycobacterium_tuberculosis.dir/
mdsTree.sh data/Mycobacterium_tuberculosis ANI 2 >& /dev/null
if ($?) exit 1
rm -r data/Mycobacterium_tuberculosis.dir/


echo ""
echo "Perfect tree ..."
makeDistTree -qc -data data/tree4 -dissim dist | grep -v '^CHRON: ' > tree4.makeDistTree
if ($?) exit 1
diff tree4.makeDistTree data/tree4.makeDistTree
if ($?) exit 1
rm tree4.makeDistTree
echo "Verbose..."
makeDistTree -qc -data data/tree4 -dissim dist  >& /dev/null
if ($?) exit 1


echo ""
echo "mdsTree: Random tree ..."
rm -r -f data/randomTree.dir/
mdsTree.sh data/randomTree dist 2 >& /dev/null
if ($?) exit 1
makeDistTree  -qc  -input_tree data/randomTree.dir/  -data data/randomTree  -dissim dist  -variance lin  -optimize  -output_tree random-output.tree > /dev/null
if ($?) exit 1
makeDistTree  -qc  -input_tree random-output.tree  -data data/randomTree  -dissim dist  -variance lin | grep -v '^CHRON: ' > randomTree.makeDistTree
if ($?) exit 1
diff randomTree.makeDistTree data/randomTree.makeDistTree
if ($?) exit 1
echo "Verbose..."
makeDistTree  -qc  -verbose 2  -input_tree random-output.tree  -data data/randomTree  -dissim dist  -variance lin  >& /dev/null
if ($?) exit 1
rm -r data/randomTree.dir/
rm randomTree.makeDistTree
rm random-output.tree


echo ""
echo "prot-identical_comm: subgraphs ..."
# Check time ??
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  \
  -optimize  \
  -delete_outliers prot-identical_comm.outliers \
  -find_hybrids prot-identical_comm.hybrids  -delete_hybrids \
  | grep -v '^CHRON: ' > prot-identical_comm.distTree
if ($?) exit 1
diff prot-identical_comm.outliers data/prot-identical_comm.outliers
if ($?) exit 1
rm prot-identical_comm.outliers
diff prot-identical_comm.hybrids data/prot-identical_comm.hybrids
if ($?) exit 1
rm prot-identical_comm.hybrids
diff prot-identical_comm.distTree data/prot-identical_comm.distTree
if ($?) exit 1
rm prot-identical_comm.distTree

echo ""
echo "prot-identical_comm: subgraphs, threads ..."
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  -optimize  -threads 3 > /dev/null
if ($?) exit 1

echo ""
echo "prot-identical_comm: subgraphs, delete ..."
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  -delete data/delete.list > /dev/null
if ($?) exit 1

echo ""
echo ""
echo "Saccharomyces hybrids ..."
cp data/Saccharomyces.dm .
if ($?) exit 1
cp /dev/null Saccharomyces.hybrid.all
while (1)
	echo ""
	makeDistTree  -threads 3  -data Saccharomyces  -dissim cons  -optimize  -subgraph_iter_max 3  -find_hybrids Saccharomyces.hybrid > /dev/null
	if ($?) exit 1
	wc -l Saccharomyces.hybrid
	if (-z Saccharomyces.hybrid) break
	cat Saccharomyces.hybrid | sort >> Saccharomyces.hybrid.all
	if ($?) exit 1
	cut -f 1 Saccharomyces.hybrid > Saccharomyces.hybrid.list
	if ($?) exit 1
	dm2subset Saccharomyces Saccharomyces.hybrid.list -exclude > Saccharomyces1.dm
	if ($?) exit 1
	rm Saccharomyces.hybrid.list
	mv Saccharomyces1.dm Saccharomyces.dm
	if ($?) exit 1
end
rm Saccharomyces.hybrid
diff Saccharomyces.hybrid.all data/Saccharomyces.hybrid.all
if ($?) exit 1
rm Saccharomyces.hybrid.all
rm Saccharomyces.dm

echo ""
echo ""
echo "prot-identical_comm: whole ..."
# Time: 7 min.
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  \
  -optimize  -whole  \
  -delete_outliers prot-identical_comm.outliers-whole \
  | grep -v '^CHRON: ' > prot-identical_comm.distTree-whole
if ($?) exit 1
diff prot-identical_comm.outliers-whole data/prot-identical_comm.outliers-whole
if ($?) exit 1
rm prot-identical_comm.outliers-whole
diff prot-identical_comm.distTree-whole data/prot-identical_comm.distTree-whole
if ($?) exit 1
rm prot-identical_comm.distTree-whole


# ??
if (0) then
  makeDistTree -data data/sample299-prot_core -dissim cons -variance lin -optimize -whole
  # was: absCriterion = 75.64
endif

