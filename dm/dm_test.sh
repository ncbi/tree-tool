#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Data Master test"
  echo "#1: seed (>=1)"
  exit 1
fi


echo ""
echo "numeric"
$THIS/numeric_test -qc go 

echo ""
echo "matrix"
$THIS/matrix_test -qc go

echo ""
echo "dataset"
$THIS/dataset_test -qc  -seed $1  go

echo ""
echo "pca: pca"
# 2 1 0
$THIS/pca  -qc  -maxTotalExpl 1  -minExpl 0  data/pca pca  
rm data/pca-pca.*

echo ""
echo "pca: Listeria_monocytogenes-qc"
cp data/Listeria_monocytogenes-qc.dm .
# -attr_pvalue 1e-5  5 1 0.07
$THIS/pca  -qc  -maxClusters 4  -mds  Listeria_monocytogenes-qc pc 
diff Listeria_monocytogenes-qc-pc.txt data/Listeria_monocytogenes-qc-pc.txt
diff Listeria_monocytogenes-qc-pc.dm data/Listeria_monocytogenes-qc-pc.dm
diff Listeria_monocytogenes-qc-pc.mds data/Listeria_monocytogenes-qc-pc.mds
rm Listeria_monocytogenes-qc*

echo ""
echo "mds: cities"
$THIS/mds  -qc  -attrType 1  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Dist  data/cities > cities.mds
diff cities.mds data/cities.mds
rm cities.mds

echo ""
echo "mds: blaLUT"
$THIS/mds  -qc  -attrType 0  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Similarity  data/blaLUT  > blaLUT.mds
diff blaLUT.mds data/blaLUT.mds
rm blaLUT.mds

echo ""
echo "mds: Enterobacteriaceae"
$THIS/mds  -qc  -attrType 2  -maxClusters 5  -attr Conservation  ../phylogeny/data/Enterobacteriaceae  > Enterobacteriaceae.mds
diff Enterobacteriaceae.mds data/Enterobacteriaceae.mds
rm Enterobacteriaceae.mds

echo ""
echo "mds: Peptostreptococcaceae"
$THIS/mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  data/Peptostreptococcaceae  > Peptostreptococcaceae.mds
diff Peptostreptococcaceae.mds data/Peptostreptococcaceae.mds
rm Peptostreptococcaceae.mds

echo ""
echo "mds: Bacteria"
$THIS/mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  data/Bacteria  > Bacteria.mds
diff Bacteria.mds data/Bacteria.mds
rm Bacteria.mds

echo ""
echo "mds: bla-A"
$THIS/mds  -qc  -attrType 0  -maxClusters 4  -class Class  -attr Similarity  data/bla-A  > bla-A.mds
diff bla-A.mds data/bla-A.mds
rm bla-A.mds

echo ""
echo "clust"
$THIS/clust  -qc  data/hmmScore 100 10 0.5 > hmmScore.clust
diff hmmScore.clust data/hmmScore.clust
rm hmmScore.clust

echo ""
echo "linreg"
$THIS/linreg_test  -qc  -seed $1
$THIS/linreg_test  -qc  -lin_dep  -seed $1

echo ""
echo "logreg"
$THIS/logreg_test  -qc  -seed $1
$THIS/logreg_test  -qc  -lin_dep  -seed $1

echo ""
echo "logreg GENOME.dm"
$THIS/logreg  -qc  data/GENOME target > logreg.out
diff logreg.out data/GENOME.logreg
rm logreg.out

echo ""
echo "Beta1"
$THIS/beta1_test  -qc  1000  -seed $1 

#echo ""
#echo "Zipf" ??
#$THIS/testZipf  -qc  1000  -seed $1

echo ""
echo "uniKernel"
$THIS/uniKernel  -qc  data/bimodal Z1_1 > bimodal.uniKernel
diff bimodal.uniKernel data/bimodal.uniKernel
rm bimodal.uniKernel

echo ""
$THIS/uniKernel  -qc  data/363068-2319168-diff JK > 363068-2319168-diff.uniKernel
diff 363068-2319168-diff.uniKernel data/363068-2319168-diff.uniKernel
rm 363068-2319168-diff.uniKernel

echo ""
echo "clust_binomial"
$THIS/clust_binomial -qc  data/Fungi-univ-stat V1 273 -outlier_pValue 1e-5 > Fungi-univ-stat.clust_binomial
diff Fungi-univ-stat.clust_binomial data/Fungi-univ-stat.clust_binomial
rm Fungi-univ-stat.clust_binomial

