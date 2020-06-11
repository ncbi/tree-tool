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
$THIS/pca  -qc  -maxTotalExpl 1  -minExpl 0  $THIS/data/pca pca  
rm $THIS/data/pca-pca.*

echo ""
echo "pca: Listeria_monocytogenes-qc"
cp $THIS/data/Listeria_monocytogenes-qc.dm .
# -attr_pvalue 1e-5  5 1 0.07
$THIS/pca  -qc  -maxClusters 4  -mds  Listeria_monocytogenes-qc pc 
diff Listeria_monocytogenes-qc-pc.txt $THIS/data/Listeria_monocytogenes-qc-pc.txt
diff Listeria_monocytogenes-qc-pc.dm $THIS/data/Listeria_monocytogenes-qc-pc.dm
diff Listeria_monocytogenes-qc-pc.mds $THIS/data/Listeria_monocytogenes-qc-pc.mds
rm Listeria_monocytogenes-qc*

echo ""
echo "mds: cities"
$THIS/mds  -qc  -attrType 1  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Dist  $THIS/data/cities > cities.mds
diff cities.mds $THIS/data/cities.mds
rm cities.mds

echo ""
echo "mds: blaLUT"
$THIS/mds  -qc  -attrType 0  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Similarity  $THIS/data/blaLUT  > blaLUT.mds
diff blaLUT.mds $THIS/data/blaLUT.mds
rm blaLUT.mds

echo ""
echo "mds: Enterobacteriaceae"
$THIS/mds  -qc  -attrType 2  -maxClusters 5  -attr Conservation  $THIS/../phylogeny/data/Enterobacteriaceae  > Enterobacteriaceae.mds
diff Enterobacteriaceae.mds $THIS/data/Enterobacteriaceae.mds
rm Enterobacteriaceae.mds

echo ""
echo "mds: Peptostreptococcaceae"
$THIS/mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  $THIS/data/Peptostreptococcaceae  > Peptostreptococcaceae.mds
diff Peptostreptococcaceae.mds $THIS/data/Peptostreptococcaceae.mds
rm Peptostreptococcaceae.mds

echo ""
echo "mds: Bacteria"
$THIS/mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  $THIS/data/Bacteria  > Bacteria.mds
diff Bacteria.mds $THIS/data/Bacteria.mds
rm Bacteria.mds

echo ""
echo "mds: bla-A"
$THIS/mds  -qc  -attrType 0  -maxClusters 4  -class Class  -attr Similarity  $THIS/data/bla-A  > bla-A.mds
diff bla-A.mds $THIS/data/bla-A.mds
rm bla-A.mds

echo ""
echo "mdsTree: Mycobacterium_tuberculosis"
rm -rf $THIS/data/Mycobacterium_tuberculosis.dir/
$THIS/mdsTree.sh $THIS/data/Mycobacterium_tuberculosis ANI 2 &> /dev/null
rm -r $THIS/data/Mycobacterium_tuberculosis.dir/

echo ""
echo "clust"
$THIS/clust  -qc  $THIS/data/hmmScore 100 10 0.5  -threshold_SDs 3 > hmmScore.clust
diff hmmScore.clust $THIS/data/hmmScore.clust
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
$THIS/logreg  -qc  $THIS/data/GENOME target > logreg.out
diff logreg.out $THIS/data/GENOME.logreg
rm logreg.out

echo ""
echo "Beta1"
$THIS/beta1_test  -qc  1000  -seed $1 

#echo ""
#echo "Zipf" ??
#$THIS/testZipf  -qc  1000  -seed $1

echo ""
echo "uniKernel"
$THIS/uniKernel  -qc  $THIS/data/bimodal Z1_1 > bimodal.uniKernel
diff bimodal.uniKernel $THIS/data/bimodal.uniKernel
rm bimodal.uniKernel

echo ""
$THIS/uniKernel  -qc  $THIS/data/363068-2319168-diff JK > 363068-2319168-diff.uniKernel
diff 363068-2319168-diff.uniKernel $THIS/data/363068-2319168-diff.uniKernel
rm 363068-2319168-diff.uniKernel

echo ""
echo "clust_binomial"
$THIS/clust_binomial -qc  $THIS/data/Fungi-univ-stat V1 273 -outlier_pValue 1e-5 > Fungi-univ-stat.clust_binomial
diff Fungi-univ-stat.clust_binomial $THIS/data/Fungi-univ-stat.clust_binomial
rm Fungi-univ-stat.clust_binomial


echo ""
echo -e ${GREEN}SUCCESS!${NOCOLOR}
