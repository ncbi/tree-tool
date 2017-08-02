#!/bin/csh -f

if ($# != 2) then
  echo "Make a distance tree"
  echo "Output: #1.tree"
  echo "#1: Input .dm file without .dm"
  echo "#2: Distance attribute in #1"
  exit 1
endif



attr2_2phylip $1 $2 $1.map > $1.phylip
if ($?) exit 1

/panfs/pan1.be-md.ncbi.nlm.nih.gov/gpipe/ThirdParty/FastME/production/arch/x86_64/bin/fastme -i $1.phylip -o $1.phylip_tree -n -s
if ($?) exit 1
rm $1.phylip
rm $1.phylip_fastme_stat.txt

~jcherry/bin/nw_rename -l $1.phylip_tree $1.map > $1.nw
if ($?) exit 1
rm $1.phylip_tree
rm $1.map

newick2tree $1.nw > $1.tree
if ($?) exit 1
rm $1.nw

