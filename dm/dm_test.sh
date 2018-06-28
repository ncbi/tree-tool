#!/bin/csh -f


if ($# != 1) then
  echo "Data Master test"
  echo "#1: seed (>=1)"
  exit 1
endif


echo ""
echo "numeric"
numeric_test -qc go 
if ($?) exit 1

echo ""
echo "matrix"
matrix_test -qc go
if ($?) exit 1

echo ""
echo "dataset"
dataset_test -qc  -seed $1  go
if ($?) exit 1

echo ""
echo "pca: pca"
# 2 1 0
pca  -qc  -maxTotalExpl 1  -minExpl 0  data/pca pca  
if ($?) exit 1
rm data/pca-pca.*

echo ""
echo "pca: Listeria_monocytogenes-qc"
cp data/Listeria_monocytogenes-qc.dm .
# -attr_pvalue 1e-5  5 1 0.07
pca  -qc  -maxClusters 4  -mds  Listeria_monocytogenes-qc pc 
if ($?) exit 1
diff Listeria_monocytogenes-qc-pc.txt data/Listeria_monocytogenes-qc-pc.txt
if ($?) exit 1
diff Listeria_monocytogenes-qc-pc.dm data/Listeria_monocytogenes-qc-pc.dm
if ($?) exit 1
diff Listeria_monocytogenes-qc-pc.mds data/Listeria_monocytogenes-qc-pc.mds
if ($?) exit 1
rm Listeria_monocytogenes-qc*

echo ""
echo "mds: cities"
mds  -qc  -attrType 1  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Dist  data/cities > cities.mds
if ($?) exit 1
diff cities.mds data/cities.mds
if ($?) exit 1
rm cities.mds

echo ""
echo "mds: blaLUT"
mds  -qc  -attrType 0  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Similarity  data/blaLUT  > blaLUT.mds
if ($?) exit 1
diff blaLUT.mds data/blaLUT.mds
if ($?) exit 1
rm blaLUT.mds

echo ""
echo "mds: Enterobacteriaceae"
mds  -qc  -attrType 2  -maxClusters 5  -attr Conservation  ../phylogeny/data/Enterobacteriaceae  > Enterobacteriaceae.mds
if ($?) exit 1
diff Enterobacteriaceae.mds data/Enterobacteriaceae.mds
if ($?) exit 1
rm Enterobacteriaceae.mds

echo ""
echo "mds: Peptostreptococcaceae"
mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  data/Peptostreptococcaceae  > Peptostreptococcaceae.mds
if ($?) exit 1
diff Peptostreptococcaceae.mds data/Peptostreptococcaceae.mds
if ($?) exit 1
rm Peptostreptococcaceae.mds

echo ""
echo "mds: Bacteria"
mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  data/Bacteria  > Bacteria.mds
if ($?) exit 1
diff Bacteria.mds data/Bacteria.mds
if ($?) exit 1
rm Bacteria.mds

echo ""
echo "mds: bla-A"
mds  -qc  -attrType 0  -maxClusters 4  -class Class  -attr Similarity  data/bla-A  > bla-A.mds
if ($?) exit 1
diff bla-A.mds data/bla-A.mds
if ($?) exit 1
rm bla-A.mds

echo ""
echo "clust"
clust  -qc  data/hmmScore 100 10 0.5 > hmmScore.clust
if ($?) exit 1
diff hmmScore.clust data/hmmScore.clust
if ($?) exit 1
rm hmmScore.clust

echo ""
echo "linreg"
linreg_test  -qc  -seed $1
if ($?) exit 1
linreg_test  -qc  -lin_dep  -seed $1
if ($?) exit 1

echo ""
echo "logreg"
logreg_test  -qc  -seed $1
if ($?) exit 1
logreg_test  -qc  -lin_dep  -seed $1
if ($?) exit 1

echo ""
echo "logreg GENOME.dm"
logreg  -qc  data/GENOME target > logreg.out
if ($?) exit 1
diff logreg.out data/GENOME.logreg
if ($?) exit 1
rm logreg.out

echo ""
echo "Beta1"
beta1_test  -qc  1000  -seed $1 
if ($?) exit 1

#echo ""
#echo "Zipf" ??
#testZipf  -qc  1000  -seed $1
#if ($?) exit 1

echo ""
echo "uniKernel"
uniKernel  -qc  data/bimodal Z1_1 > bimodal.unikernel
if ($?) exit 1
diff bimodal.unikernel data/bimodal.unikernel
if ($?) exit 1
rm bimodal.unikernel

echo ""
uniKernel  -qc  data/363068-2319168-diff JK > 363068-2319168-diff.uniKernel
if ($?) exit 1
diff 363068-2319168-diff.uniKernel data/363068-2319168-diff.uniKernel
if ($?) exit 1
rm 363068-2319168-diff.uniKernel

echo ""
echo "clust_binomial"
clust_binomial data/Fungi-univ-stat V1 273 -outlier_pValue 1e-5 > Fungi-univ-stat.clust_binomial
if ($?) exit 1
diff Fungi-univ-stat.clust_binomial data/Fungi-univ-stat.clust_binomial
if ($?) exit 1
rm Fungi-univ-stat.clust_binomial

