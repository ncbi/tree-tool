include ../MakeRules
CPPFLAGS=$(OptCPPFLAGS) 
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
  attr2_2attr1 \
  attr2_2pairs \
  attr2_2power \
  attr2_multiply \
  attr2_nan \
  attr2_set \
  attr2log \
  attr2power \
  attrs \
  beta1_test \
  binomialCdf \
  char2attr \
  clust \
  clust_binomial \
  count \
  dataset_test \
  distTriangle \
  dm2objs \
  dm2space \
  dm2subset \
  dm_merge \
  linreg \
  linreg_test \
  logreg \
  logreg_test \
  matrix_test \
	mds \
	normalCdf \
	numeric_test \
	obj_attr2dm \
	pairs2dm \
	pca \
	pca2clusters \
	positiveAverage \
	printDataset \
	sim2dist \
  threshold \
	uniKernel \
	zipf_test
	

attr2_2attr1.o:  $(DM_HPP) 
attr2_2attr1OBJS=attr2_2attr1.o $(DM_OBJ) 
attr2_2attr1:	$(attr2_2attr1OBJS)
	$(CXX) -o $@ $(attr2_2attr1OBJS) $(LIBS)
	$(ECHO)

attr2_2pairs.o:  $(DM_HPP) 
attr2_2pairsOBJS=attr2_2pairs.o $(DM_OBJ) 
attr2_2pairs:	$(attr2_2pairsOBJS)
	$(CXX) -o $@ $(attr2_2pairsOBJS) $(LIBS)
	$(ECHO)

attr2_2power.o:  $(DM_HPP) 
attr2_2powerOBJS=attr2_2power.o $(DM_OBJ) 
attr2_2power:	$(attr2_2powerOBJS)
	$(CXX) -o $@ $(attr2_2powerOBJS) $(LIBS)
	$(ECHO)

attr2_multiply.o:  $(DM_HPP) 
attr2_multiplyOBJS=attr2_multiply.o $(DM_OBJ) 
attr2_multiply:	$(attr2_multiplyOBJS)
	$(CXX) -o $@ $(attr2_multiplyOBJS) $(LIBS)
	$(ECHO)

attr2_nan.o:  $(DM_HPP) 
attr2_nanOBJS=attr2_nan.o $(DM_OBJ) 
attr2_nan:	$(attr2_nanOBJS)
	$(CXX) -o $@ $(attr2_nanOBJS) $(LIBS)
	$(ECHO)

attr2_set.o:  $(DM_HPP) 
attr2_setOBJS=attr2_set.o $(DM_OBJ) 
attr2_set:	$(attr2_setOBJS)
	$(CXX) -o $@ $(attr2_setOBJS) $(LIBS)
	$(ECHO)

attr2log.o:  $(DM_HPP) 
attr2logOBJS=attr2log.o $(DM_OBJ) 
attr2log:	$(attr2logOBJS)
	$(CXX) -o $@ $(attr2logOBJS) $(LIBS)
	$(ECHO)

attr2power.o:  $(DM_HPP) 
attr2powerOBJS=attr2power.o $(DM_OBJ) 
attr2power:	$(attr2powerOBJS)
	$(CXX) -o $@ $(attr2powerOBJS) $(LIBS)
	$(ECHO)

attrs.o:  $(DM_HPP) 
attrsOBJS=attrs.o $(DM_OBJ) 
attrs:	$(attrsOBJS)
	$(CXX) -o $@ $(attrsOBJS) $(LIBS)
	$(ECHO)

beta1_test.o:  $(DM_HPP) 
beta1_testOBJS=beta1_test.o $(DM_OBJ)
beta1_test:	$(beta1_testOBJS)
	$(CXX) -o $@ $(beta1_testOBJS) $(LIBS)
	$(ECHO)

binomialCdf.o:  $(DM_HPP) 
binomialCdfOBJS=binomialCdf.o $(DM_OBJ)
binomialCdf:	$(binomialCdfOBJS)
	$(CXX) -o $@ $(binomialCdfOBJS) $(LIBS)
	$(ECHO)

char2attr.o:  $(DM_HPP) 
char2attrOBJS=char2attr.o $(DM_OBJ) 
char2attr:	$(char2attrOBJS)
	$(CXX) -o $@ $(char2attrOBJS) $(LIBS)
	$(ECHO)

circuit.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp $(DM_DIR)/matrix.hpp 
circuitOBJS=circuit.o $(NUMERIC_OBJ) $(DM_DIR)/matrix.o
circuit:	$(circuitOBJS)
	$(CXX) -o $@ $(circuitOBJS) $(LIBS)  
	$(ECHO)

clust.o:	$(DM_HPP) 
clustOBJS=clust.o $(DM_OBJ) 
clust:	$(clustOBJS)
	$(CXX) -o $@ $(clustOBJS) $(LIBS)
	$(ECHO)

clust_binomial.o:  $(DM_HPP) 
clust_binomialOBJS=clust_binomial.o $(DM_OBJ) 
clust_binomial:	$(clust_binomialOBJS)
	$(CXX) -o $@ $(clust_binomialOBJS) $(LIBS)
	$(ECHO)

# ??
convex.o:  $(PREDICT_HPP) 
convexOBJS=convex.o $(PREDICT_OBJ) 
convex:	$(convexOBJS)
	$(CXX) -o $@ $(convexOBJS) $(LIBS)
	$(ECHO)

count.o:	$(DM_DIR)/numeric.hpp $(COMMON_HPP) 
countOBJS=count.o $(NUMERIC_OBJ) 
count:	$(countOBJS)
	$(CXX) -o $@ $(countOBJS) $(LIBS)
	$(ECHO)

dataset_test.o:  $(DM_HPP) 
dataset_testOBJS=dataset_test.o $(DM_OBJ)
dataset_test:	$(dataset_testOBJS)
	$(CXX) -o $@ $(dataset_testOBJS) $(LIBS)
	$(ECHO)

distTriangle.o:  $(DM_HPP) 
distTriangleOBJS=distTriangle.o $(DM_OBJ) 
distTriangle:	$(distTriangleOBJS)
	$(CXX) -o $@ $(distTriangleOBJS) $(LIBS)
	$(ECHO)

dm2objs.o:  $(DM_HPP) 
dm2objsOBJS=dm2objs.o $(DM_OBJ) 
dm2objs:	$(dm2objsOBJS)
	$(CXX) -o $@ $(dm2objsOBJS) $(LIBS)
	$(ECHO)

dm2space.o:  $(DM_HPP) 
dm2spaceOBJS=dm2space.o $(DM_OBJ) 
dm2space:	$(dm2spaceOBJS)
	$(CXX) -o $@ $(dm2spaceOBJS) $(LIBS)
	$(ECHO)

dm2subset.o:  $(DM_HPP) 
dm2subsetOBJS=dm2subset.o $(DM_OBJ) 
dm2subset:	$(dm2subsetOBJS)
	$(CXX) -o $@ $(dm2subsetOBJS) $(LIBS)
	$(ECHO)

dm_merge.o:  $(DM_HPP) 
dm_mergeOBJS=dm_merge.o $(DM_OBJ) 
dm_merge:	$(dm_mergeOBJS)
	$(CXX) -o $@ $(dm_mergeOBJS) $(LIBS)
	$(ECHO)

dnf.o:  $(PREDICT_HPP) 
dnfOBJS=dnf.o $(PREDICT_OBJ)
dnf:	$(dnfOBJS)
	$(CXX) -o $@ $(dnfOBJS) $(LIBS)
	$(ECHO)

linreg.o:  $(PREDICT_HPP) 
linregOBJS=linreg.o $(PREDICT_OBJ)
linreg:	$(linregOBJS)
	$(CXX) -o $@ $(linregOBJS) $(LIBS)
	$(ECHO)

linreg_test.o:  $(PREDICT_HPP) 
linreg_testOBJS=linreg_test.o $(PREDICT_OBJ)
linreg_test:	$(linreg_testOBJS)
	$(CXX) -o $@ $(linreg_testOBJS) $(LIBS)
	$(ECHO)

logreg.o:  $(PREDICT_HPP) 
logregOBJS=logreg.o $(PREDICT_OBJ)
logreg:	$(logregOBJS)
	$(CXX) -o $@ $(logregOBJS) $(LIBS)
	$(ECHO)

logreg_test.o:  $(PREDICT_HPP) 
logreg_testOBJS=logreg_test.o $(PREDICT_OBJ)
logreg_test:	$(logreg_testOBJS)
	$(CXX) -o $@ $(logreg_testOBJS) $(LIBS)
	$(ECHO)

matrix_test.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp $(DM_DIR)/matrix.hpp 
matrix_testOBJS=matrix_test.o $(NUMERIC_OBJ) $(DM_DIR)/matrix.o
matrix_test:	$(matrix_testOBJS)
	$(CXX) -o $@ $(matrix_testOBJS) $(LIBS)  
	$(ECHO)

mds.o:  $(PREDICT_HPP) 
mdsOBJS=mds.o $(PREDICT_OBJ)
mds:	$(mdsOBJS)
	$(CXX) -o $@ $(mdsOBJS) $(LIBS)
	$(ECHO)

normalCdf.o:  $(DM_HPP) 
normalCdfOBJS=normalCdf.o $(DM_OBJ)
normalCdf:	$(normalCdfOBJS)
	$(CXX) -o $@ $(normalCdfOBJS) $(LIBS)
	$(ECHO)

numeric_test.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp 
numeric_testOBJS=numeric_test.o $(NUMERIC_OBJ)
numeric_test:	$(numeric_testOBJS)
	$(CXX) -o $@ $(numeric_testOBJS) $(LIBS)
	$(ECHO)

obj_attr2dm.o:  $(DM_HPP) 
obj_attr2dmOBJS=obj_attr2dm.o $(DM_OBJ) 
obj_attr2dm:	$(obj_attr2dmOBJS)
	$(CXX) -o $@ $(obj_attr2dmOBJS) $(LIBS)
	$(ECHO)

pairs2dm.o:  $(DM_HPP) 
pairs2dmOBJS=pairs2dm.o $(DM_OBJ) 
pairs2dm:	$(pairs2dmOBJS)
	$(CXX) -o $@ $(pairs2dmOBJS) $(LIBS)
	$(ECHO)

pca.o:  $(PREDICT_HPP) 
pcaOBJS=pca.o $(PREDICT_OBJ)
pca:	$(pcaOBJS)
	$(CXX) -o $@ $(pcaOBJS) $(LIBS)
	$(ECHO)

pca2clusters.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp 
pca2clustersOBJS=pca2clusters.o $(NUMERIC_OBJ)
pca2clusters:	$(pca2clustersOBJS)
	$(CXX) -o $@ $(pca2clustersOBJS) $(LIBS)
	$(ECHO)

positiveAverage.o:  $(DM_HPP) 
positiveAverageOBJS=positiveAverage.o $(DM_OBJ)
positiveAverage:	$(positiveAverageOBJS)
	$(CXX) -o $@ $(positiveAverageOBJS) $(LIBS)
	$(ECHO)

printDataset.o:  $(DM_HPP) 
printDatasetOBJS=printDataset.o $(DM_OBJ) 
printDataset:	$(printDatasetOBJS)
	$(CXX) -o $@ $(printDatasetOBJS) $(LIBS)
	$(ECHO)

sim2dist.o:  $(DM_HPP) 
sim2distOBJS=sim2dist.o $(DM_OBJ) 
sim2dist:	$(sim2distOBJS)
	$(CXX) -o $@ $(sim2distOBJS) $(LIBS)
	$(ECHO)

threshold.o:  $(DM_HPP) 
thresholdOBJS=threshold.o $(DM_OBJ) 
threshold:	$(thresholdOBJS)
	$(CXX) -o $@ $(thresholdOBJS) $(LIBS)
	$(ECHO)

uniKernel.o:  $(DM_HPP) 
uniKernelOBJS=uniKernel.o $(DM_OBJ) 
uniKernel:	$(uniKernelOBJS)
	$(CXX) -o $@ $(uniKernelOBJS) $(LIBS)
	$(ECHO)

zipf_test.o:  $(DM_HPP) 
zipf_testOBJS=zipf_test.o $(DM_OBJ) 
zipf_test:	$(zipf_testOBJS)
	$(CXX) -o $@ $(zipf_testOBJS) $(LIBS)
	$(ECHO)

