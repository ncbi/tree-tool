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

makeDistTree -qc -input_tree data/Enterobacteriaceae.dir/ -data data/Enterobacteriaceae -dissim Conservation -topology -output_tree Enterobacteriaceae.tree > /dev/null
if ($?) exit 1
rm -r data/Enterobacteriaceae.dir/

printDistTree -qc -data data/Enterobacteriaceae -dissim Conservation Enterobacteriaceae.tree > Enterobacteriaceae.nw
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
makeDistTree  -qc  -input_tree data/randomTree.dir/  -data data/randomTree  -dissim dist  -variance lin  -topology  -output_tree random-output.tree >& /dev/null
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
echo "prot-identical_comm ..."
# Check time ??
makeDistTree  -qc  -data data/prot-identical_comm  -dissim cons  -topology  -remove_outliers prot-identical_comm.outliers | grep -v '^CHRON: ' > prot-identical_comm.distTree
if ($?) exit 1
diff prot-identical_comm.outliers data/prot-identical_comm.outliers
if ($?) exit 1
rm prot-identical_comm.outliers
diff prot-identical_comm.distTree data/prot-identical_comm.distTree
if ($?) exit 1
rm prot-identical_comm.distTree




