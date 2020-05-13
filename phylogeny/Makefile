include ../MakeRules
CPPFLAGS=$(OptCPPFLAGS) 
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
  attr2_2paup \
  attr2_2phylip \
  combine_dissims \
  compareTrees \
  distTree_new \
  feature2dissim \
  feature2gain_loss \
  hash2dissim \
  hash_request2dissim \
  loci_request2dissim \
	makeDistTree \
	makeDistTree-dev \
	makeFeatureTree \
	mergeSimilarity \
	mlst2dissim \
	mutation2phen \
	newick2tree \
	printDistTree \
	randomDistTree \
  sampleReport \
	statDistTree \
	testDistTree \
	tree2clusters \
	tree2genogroup
	

attr2_2paup.o:  $(DM_HPP) $(PHYL_DIR)/version.inc
attr2_2paupOBJS=attr2_2paup.o $(DM_OBJ) 
attr2_2paup:	$(attr2_2paupOBJS)
	$(CXX) -o $@ $(attr2_2paupOBJS) $(LIBS)
	$(ECHO)

attr2_2phylip.o:  $(DM_HPP) $(PHYL_DIR)/version.inc
attr2_2phylipOBJS=attr2_2phylip.o $(DM_OBJ) 
attr2_2phylip:	$(attr2_2phylipOBJS)
	$(CXX) -o $@ $(attr2_2phylipOBJS) $(LIBS)
	$(ECHO)

attr2_multiply.o:  $(DM_HPP) $(PHYL_DIR)/version.inc
attr2_multiplyOBJS=attr2_multiply.o $(DM_OBJ) 
attr2_multiply:	$(attr2_multiplyOBJS)
	$(CXX) -o $@ $(attr2_multiplyOBJS) $(LIBS)
	$(ECHO)

combine_dissims.o:  $(NUMERIC_HPP) $(PHYL_DIR)/version.inc
combine_dissimsOBJS=combine_dissims.o $(CPP_DIR)/common.o $(DM_DIR)/numeric.o
combine_dissims:	$(combine_dissimsOBJS)
	$(CXX) -o $@ $(combine_dissimsOBJS) $(LIBS)
	$(ECHO)

compareTrees.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/featureTree.hpp $(PHYL_DIR)/version.inc
compareTreesOBJS=compareTrees.o $(DISTTREE_OBJ) $(PHYL_DIR)/featureTree.o
compareTrees:	$(compareTreesOBJS)
	$(CXX) -o $@ $(compareTreesOBJS) $(LIBS)
	$(ECHO)
	
distTree_new.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
distTree_newOBJS=distTree_new.o $(DISTTREE_OBJ)
distTree_new:	$(distTree_newOBJS)
	$(CXX) -o $@ $(distTree_newOBJS) $(LIBS)
	$(ECHO)

feature2dissim.o:  $(DM_HPP) $(PHYL_DIR)/version.inc
feature2dissimOBJS=feature2dissim.o $(DM_OBJ) 
feature2dissim:	$(feature2dissimOBJS)
	$(CXX) -o $@ $(feature2dissimOBJS) $(LIBS)
	$(ECHO)

feature2gain_loss.o:  $(COMMON_HPP) $(PHYL_DIR)/featureTree.hpp $(PHYL_DIR)/version.inc
feature2gain_lossOBJS=feature2gain_loss.o $(FEATURE_TREE_OBJ)
feature2gain_loss:	$(feature2gain_lossOBJS)
	$(CXX) -o $@ $(feature2gain_lossOBJS) $(LIBS)
	$(ECHO)

hash2dissim.o:  $(DM_HPP) $(PHYL_DIR)/evolution.hpp $(PHYL_DIR)/version.inc
hash2dissimOBJS=hash2dissim.o $(DM_OBJ) $(PHYL_DIR)/evolution.o
hash2dissim:	$(hash2dissimOBJS)
	$(CXX) -o $@ $(hash2dissimOBJS) $(LIBS)
	$(ECHO)

hash_request2dissim.o:  $(NUMERIC_HPP) $(PHYL_DIR)/evolution.hpp $(PHYL_DIR)/version.inc
hash_request2dissimOBJS=hash_request2dissim.o $(DM_OBJ) $(PHYL_DIR)/evolution.o
hash_request2dissim:	$(hash_request2dissimOBJS)
	$(CXX) -o $@ $(hash_request2dissimOBJS) $(LIBS)
	$(ECHO)

loci_request2dissim.o:  $(COMMON_HPP) $(PHYL_DIR)/version.inc
loci_request2dissimOBJS=loci_request2dissim.o $(CPP_DIR)/common.o
loci_request2dissim:	$(loci_request2dissimOBJS)
	$(CXX) -o $@ $(loci_request2dissimOBJS) $(LIBS)
	$(ECHO)

makeDistTree.o:  $(COMMON_HPP) $(CPP_DIR)/graph.hpp $(DM_DIR)/dataset.hpp $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
makeDistTreeOBJS=makeDistTree.o $(DISTTREE_OBJ)
makeDistTree:	$(makeDistTreeOBJS)
	$(CXX) -o $@ $(makeDistTreeOBJS) $(LIBS)
	$(ECHO)
makeDistTree-dev:	$(makeDistTreeOBJS)
	$(CXX) -o $@ $(makeDistTreeOBJS) $(LIBS)
	$(ECHO)

makeFeatureTree.o:  $(NUMERIC_HPP) $(CPP_DIR)/graph.hpp $(PHYL_DIR)/featureTree.hpp $(PHYL_DIR)/version.inc
makeFeatureTreeOBJS=makeFeatureTree.o $(FEATURE_TREE_OBJ)
makeFeatureTree:	$(makeFeatureTreeOBJS)
	$(CXX) -o $@ $(makeFeatureTreeOBJS) $(LIBS)
	$(ECHO)

mergeSimilarity.o:  $(DM_HPP) $(PHYL_DIR)/version.inc
mergeSimilarityOBJS=mergeSimilarity.o $(DM_OBJ) 
mergeSimilarity:	$(mergeSimilarityOBJS)
	$(CXX) -o $@ $(mergeSimilarityOBJS) $(LIBS)
	$(ECHO)

mlst2dissim.o:  $(DM_HPP) $(PHYL_DIR)/version.inc
mlst2dissimOBJS=mlst2dissim.o $(DM_OBJ) 
mlst2dissim:	$(mlst2dissimOBJS)
	$(CXX) -o $@ $(mlst2dissimOBJS) $(LIBS)
	$(ECHO)

mutation2phen.o:  $(COMMON_HPP) $(PHYL_DIR)/version.inc
mutation2phenOBJS=mutation2phen.o $(CPP_DIR)/common.o
mutation2phen:	$(mutation2phenOBJS)
	$(CXX) -o $@ $(mutation2phenOBJS) $(LIBS)
	$(ECHO)

newick2tree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
newick2treeOBJS=newick2tree.o $(DISTTREE_OBJ)
newick2tree:	$(newick2treeOBJS)
	$(CXX) -o $@ $(newick2treeOBJS) $(LIBS)
	$(ECHO)

printDistTree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
printDistTreeOBJS=printDistTree.o $(DISTTREE_OBJ)
printDistTree:	$(printDistTreeOBJS)
	$(CXX) -o $@ $(printDistTreeOBJS) $(LIBS)
	$(ECHO)

randomDistTree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
randomDistTreeOBJS=randomDistTree.o $(DISTTREE_OBJ)
randomDistTree:	$(randomDistTreeOBJS)
	$(CXX) -o $@ $(randomDistTreeOBJS) $(LIBS)
	$(ECHO)

sampleReport.o:  $(COMMON_HPP) $(PHYL_DIR)/version.inc
# $(PHYL_DIR)/distTree.hpp
sampleReportOBJS=sampleReport.o $(CPP_DIR)/common.o
# $(DISTTREE_OBJ)
sampleReport:	$(sampleReportOBJS)
	$(CXX) -o $@ $(sampleReportOBJS) $(LIBS)
	$(ECHO)

statDistTree.o:  $(COMMON_HPP) $(DM_DIR)/dataset.hpp $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
statDistTreeOBJS=statDistTree.o $(DISTTREE_OBJ)
statDistTree:	$(statDistTreeOBJS)
	$(CXX) -o $@ $(statDistTreeOBJS) $(LIBS)
	$(ECHO)

testDistTree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
testDistTreeOBJS=testDistTree.o $(DISTTREE_OBJ)
testDistTree:	$(testDistTreeOBJS)
	$(CXX) -o $@ $(testDistTreeOBJS) $(LIBS)
	$(ECHO)

tree2clusters.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
tree2clustersOBJS=tree2clusters.o $(DISTTREE_OBJ)
tree2clusters:	$(tree2clustersOBJS)
	$(CXX) -o $@ $(tree2clustersOBJS) $(LIBS)
	$(ECHO)

tree2genogroup.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/version.inc
tree2genogroupOBJS=tree2genogroup.o $(DISTTREE_OBJ)
tree2genogroup:	$(tree2genogroupOBJS)
	$(CXX) -o $@ $(tree2genogroupOBJS) $(LIBS)
	$(ECHO)

