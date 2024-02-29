#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Data Master test"
  echo "#1: seed (>=1)"
  exit 1
fi
SEED=$1



TMP=`mktemp`


cd $THIS


section "numeric"
$THIS/numeric_test -qc go 

section "matrix"
$THIS/matrix_test -qc go

section "dataset"
$THIS/dataset_test -qc  -seed $SEED  go


super_section "PCA"
section "pca"
# 2 1 0
$THIS/pca  -qc  -maxTotalExpl 1  -minExpl 0  $THIS/data/pca pca  
rm $THIS/data/pca-pca.*

section "Listeria_monocytogenes-qc"
cp $THIS/data/Listeria_monocytogenes-qc.dm .
# -attr_pvalue 1e-5  5 1 0.07
$THIS/pca  -qc  -maxClusters 4  -mds  Listeria_monocytogenes-qc pc 
diff Listeria_monocytogenes-qc-pc.txt $THIS/data/Listeria_monocytogenes-qc-pc.txt
diff -b Listeria_monocytogenes-qc-pc.dm $THIS/data/Listeria_monocytogenes-qc-pc.dm
diff -b Listeria_monocytogenes-qc-pc.mds $THIS/data/Listeria_monocytogenes-qc-pc.mds
rm Listeria_monocytogenes-qc.dm
rm Listeria_monocytogenes-qc-pc.*


super_section "MDS"
section "cities"
$THIS/mds  -qc  -attrType 1  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Dist  $THIS/data/cities > cities.mds
diff cities.mds $THIS/data/cities.mds
rm cities.mds

section "blaLUT"
$THIS/mds  -qc  -attrType 0  -maxAttr 2  -maxTotalExpl 1  -minExpl 0  -attr Similarity  $THIS/data/blaLUT  > blaLUT.mds
diff blaLUT.mds $THIS/data/blaLUT.mds
rm blaLUT.mds

section "Enterobacteriaceae"
gunzip -c $THIS/../phylogeny/data/Enterobacteriaceae.dm.gz > $TMP.dm
$THIS/mds  -qc  -attrType 2  -maxClusters 5  -attr Conservation  $TMP  > Enterobacteriaceae.mds
diff Enterobacteriaceae.mds $THIS/data/Enterobacteriaceae.mds
rm Enterobacteriaceae.mds

section "Peptostreptococcaceae"
$THIS/mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  $THIS/data/Peptostreptococcaceae  > Peptostreptococcaceae.mds
diff Peptostreptococcaceae.mds $THIS/data/Peptostreptococcaceae.mds
rm Peptostreptococcaceae.mds

section "Bacteria"
$THIS/mds  -qc  -attrType 2  -maxClusters 4  -attr Conservation  $THIS/data/Bacteria  > Bacteria.mds
diff Bacteria.mds $THIS/data/Bacteria.mds
rm Bacteria.mds

section "bla-A"
$THIS/mds  -qc  -attrType 0  -maxClusters 4  -class Class  -attr Similarity  $THIS/data/bla-A  > bla-A.mds
diff bla-A.mds $THIS/data/bla-A.mds
rm bla-A.mds

section "mdsTree: Mycobacterium_tuberculosis"
rm -rf $THIS/data/Mycobacterium_tuberculosis.dir/
$THIS/mdsTree.sh $THIS/data/Mycobacterium_tuberculosis ANI 2 &> /dev/null
rm -r $THIS/data/Mycobacterium_tuberculosis.dir/


super_section "Clustering"
section "clust"
$THIS/clust  -qc  $THIS/data/hmmScore 100 10 0.5  -threshold_SDs 3 > hmmScore.clust
diff hmmScore.clust $THIS/data/hmmScore.clust
rm hmmScore.clust

section "clust_binomial"
$THIS/clust_binomial -qc  $THIS/data/Fungi-univ-stat V1 273 -outlier_pValue 1e-5 > Fungi-univ-stat.clust_binomial
diff Fungi-univ-stat.clust_binomial $THIS/data/Fungi-univ-stat.clust_binomial
rm Fungi-univ-stat.clust_binomial


super_section "Prediction"
section "linreg"
$THIS/linreg_test  -qc  -seed $SEED
$THIS/linreg_test  -qc  -lin_dep  -seed $SEED

section "logreg"
$THIS/logreg_test  -qc  -seed $SEED
$THIS/logreg_test  -qc  -lin_dep  -seed $SEED

section "logreg GENOME.dm"
$THIS/logreg  -qc  $THIS/data/GENOME target > logreg.out
diff logreg.out $THIS/data/GENOME.logreg
rm logreg.out


super_section "Distributions"
section "Beta1"
$THIS/beta1_test  -qc  1000  -seed $SEED 

#section "Zipf" ??
#$THIS/testZipf  -qc  1000  -seed $SEED

section "uniKernel"
$THIS/uniKernel  -qc  $THIS/data/bimodal Z1_1 > bimodal.uniKernel
diff bimodal.uniKernel $THIS/data/bimodal.uniKernel
rm bimodal.uniKernel

echo ""
$THIS/uniKernel  -qc  $THIS/data/363068-2319168-diff JK > 363068-2319168-diff.uniKernel
diff 363068-2319168-diff.uniKernel $THIS/data/363068-2319168-diff.uniKernel
rm 363068-2319168-diff.uniKernel


echo ""
success


rm $TMP*
