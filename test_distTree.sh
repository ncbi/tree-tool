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

makeDistTree data/Enterobacteriaceae.dir/ -data data/Enterobacteriaceae -dissim Conservation -topology -output_tree Enterobacteriaceae.tree > /dev/null
if ($?) exit 1
rm -r data/Enterobacteriaceae.dir/

printDistTree -data data/Enterobacteriaceae -dissim Conservation Enterobacteriaceae.tree > Enterobacteriaceae.nw
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
makeDistTree "" -data data/tree4 -dissim dist | grep -v '^CHRON: ' > tree4.makeDistTree
if ($?) exit 1
diff tree4.makeDistTree data/tree4.makeDistTree
if ($?) exit 1
rm tree4.makeDistTree
echo "Verbose..."
makeDistTree "" -data data/tree4 -dissim dist -verbose 2 >& /dev/null
if ($?) exit 1


echo ""
echo "mdsTree: Random tree ..."
rm -r -f data/randomTree.dir/
mdsTree.sh data/randomTree dist 2 >& /dev/null
if ($?) exit 1
makeDistTree data/randomTree.dir/ -data data/randomTree -dissim dist -variance lin -topology -output_tree random-output.tree >& /dev/null
if ($?) exit 1
makeDistTree random-output.tree -data data/randomTree -dissim dist -variance lin | grep -v '^CHRON: ' > randomTree.makeDistTree
if ($?) exit 1
diff randomTree.makeDistTree data/randomTree.makeDistTree
if ($?) exit 1
echo "Verbose..."
makeDistTree random-output.tree -data data/randomTree -dissim dist -variance lin -verbose 2  >& /dev/null
if ($?) exit 1
rm -r data/randomTree.dir/
rm randomTree.makeDistTree
rm random-output.tree


echo ""
echo "prot-identical_comm ..."
# Check time ??
distTree.sh data/prot-identical_comm cons | grep -v '^CHRON: ' > prot-identical_comm.distTree
if ($?) exit 1
diff prot-identical_comm.distTree data/prot-identical_comm.distTree
if ($?) exit 1


