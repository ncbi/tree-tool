#!/bin/csh -f

if ($# != 3) then
  echo "Hierarchical clustering by MDS. Output directory: #1/"
  echo "#1: Input .dm file without .dm"
  echo "#2: Pairwise measurement attribute in #1"
  echo "#3: mds -attrType"
  exit 1
endif


set EXEC = ~brovervv/code


set DEBUG = 0  # or 1 


# TmpFNam, Verbose
if ($DEBUG) then
  set TmpFNam = tmp
  set Verbose = 2
else
  set TmpFNam = `mktemp` 
  set Verbose = 0
endif

set exitcode = 1


mkdir $1.dir
if ($?) exit 1

mkdir $TmpFNam.tmp
if ($?) exit 1

cp $1.dm $TmpFNam.tmp/
if ($?) exit 1

echo "" > $TmpFNam.mds

while (1)
  echo ""
  set F = `ls $TmpFNam.tmp | head -1`
  if ($F == "")  break
  set F = `basename $F .dm`
  echo $F ...
  set N = `grep -v '^#' $TmpFNam.tmp/$F.dm | head -1`
  echo "# objects = " $N[2]
  if ($N[2] <= 5) then  # PAR
    mv $TmpFNam.tmp/$F.dm $1.dir/
    if ($?) exit 1
  else
    echo "" >> $TmpFNam.mds
    echo "" >> $TmpFNam.mds
    set L_old = `ls $TmpFNam.tmp/ | wc -l`
    $EXEC/cpp/dm/mds -verbose $Verbose  -attrType $3  -maxClusters 6  -clusteringDir $TmpFNam.tmp  -attr $2  $TmpFNam.tmp/$F >> $TmpFNam.mds
    if ($?) exit 1
    set L_new = `ls $TmpFNam.tmp/ | wc -l`
    if ($L_old == $L_new) then
      mv $TmpFNam.tmp/$F.dm $1.dir/
      if ($?) exit 1
    else
      rm $TmpFNam.tmp/$F.dm
      if ($?) exit 1
    endif
  endif
end


set exitcode = 0
quit:
if (! $DEBUG)  rm -fr $TmpFNam*   
exit $exitcode
