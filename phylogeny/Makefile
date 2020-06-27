include ../MakeRules
CPPFLAGS=$(OptCPPFLAGS) 
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
	asnt2tree \
  attr2_2paup \
  attr2_2phylip \
  compareTrees \
  distTree_new \
  dm2phen \
  feature2gain_loss \
	makeDistTree \
	makeDistTree-dev \
	makeFeatureTree \
	mergeSimilarity \
	mutation2phen \
	newick2tree \
	printDistTree \
	randomDistTree \
  sampleReport \
	statDistTree \
	testDistTree \
	tree2clusters \
	tree2genogroup
	

asnt2tree.o:  $(COMMON_HPP) $(CPP_DIR)/asn.hpp $(DM_DIR)/numeric.hpp $(PHYL_DIR)/distTree.hpp 
asnt2treeOBJS=asnt2tree.o $(DISTTREE_OBJ) $(CPP_DIR)/asn.o
asnt2tree:	$(asnt2treeOBJS)
	$(CXX) -o $@ $(asnt2treeOBJS) $(LIBS) 
	$(ECHO)

attr2_2paup.o:  $(DM_HPP) 
attr2_2paupOBJS=attr2_2paup.o $(DM_OBJ) 
attr2_2paup:	$(attr2_2paupOBJS)
	$(CXX) -o $@ $(attr2_2paupOBJS) $(LIBS)
	$(ECHO)

attr2_2phylip.o:  $(DM_HPP) 
attr2_2phylipOBJS=attr2_2phylip.o $(DM_OBJ) 
attr2_2phylip:	$(attr2_2phylipOBJS)
	$(CXX) -o $@ $(attr2_2phylipOBJS) $(LIBS)
	$(ECHO)

dm2phen.o:  $(DM_HPP) 
dm2phenOBJS=dm2phen.o $(DM_OBJ) 
dm2phen:	$(dm2phenOBJS)
	$(CXX) -o $@ $(dm2phenOBJS) $(LIBS)
	$(ECHO)

compareTrees.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/featureTree.hpp 
compareTreesOBJS=compareTrees.o $(DISTTREE_OBJ) $(PHYL_DIR)/featureTree.o
compareTrees:	$(compareTreesOBJS)
	$(CXX) -o $@ $(compareTreesOBJS) $(LIBS)
	$(ECHO)
	
distTree_new.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp 
distTree_newOBJS=distTree_new.o $(DISTTREE_OBJ)
distTree_new:	$(distTree_newOBJS)
	$(CXX) -o $@ $(distTree_newOBJS) $(LIBS)
	$(ECHO)

feature2gain_loss.o:  $(COMMON_HPP) $(PHYL_DIR)/featureTree.hpp 
feature2gain_lossOBJS=feature2gain_loss.o $(FEATURE_TREE_OBJ)
feature2gain_loss:	$(feature2gain_lossOBJS)
	$(CXX) -o $@ $(feature2gain_lossOBJS) $(LIBS)
	$(ECHO)

makeDistTree.o:  $(COMMON_HPP) $(CPP_DIR)/graph.hpp $(DM_DIR)/dataset.hpp $(PHYL_DIR)/distTree.hpp 
makeDistTreeOBJS=makeDistTree.o $(DISTTREE_OBJ)
makeDistTree:	$(makeDistTreeOBJS)
	$(CXX) -o $@ $(makeDistTreeOBJS) $(LIBS)
	$(ECHO)
makeDistTree-dev:	$(makeDistTreeOBJS)
	$(CXX) -o $@ $(makeDistTreeOBJS) $(LIBS)
	$(ECHO)

makeFeatureTree.o:  $(NUMERIC_HPP) $(CPP_DIR)/graph.hpp $(PHYL_DIR)/featureTree.hpp 
makeFeatureTreeOBJS=makeFeatureTree.o $(FEATURE_TREE_OBJ)
makeFeatureTree:	$(makeFeatureTreeOBJS)
	$(CXX) -o $@ $(makeFeatureTreeOBJS) $(LIBS)
	$(ECHO)

mergeSimilarity.o:  $(DM_HPP) 
mergeSimilarityOBJS=mergeSimilarity.o $(DM_OBJ) 
mergeSimilarity:	$(mergeSimilarityOBJS)
	$(CXX) -o $@ $(mergeSimilarityOBJS) $(LIBS)
	$(ECHO)

mutation2phen.o:  $(COMMON_HPP) 
mutation2phenOBJS=mutation2phen.o $(CPP_DIR)/common.o
mutation2phen:	$(mutation2phenOBJS)
	$(CXX) -o $@ $(mutation2phenOBJS) $(LIBS)
	$(ECHO)

newick2tree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp 
newick2treeOBJS=newick2tree.o $(DISTTREE_OBJ)
newick2tree:	$(newick2treeOBJS)
	$(CXX) -o $@ $(newick2treeOBJS) $(LIBS)
	$(ECHO)

printDistTree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp 
printDistTreeOBJS=printDistTree.o $(DISTTREE_OBJ)
printDistTree:	$(printDistTreeOBJS)
	$(CXX) -o $@ $(printDistTreeOBJS) $(LIBS)
	$(ECHO)

randomDistTree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp 
randomDistTreeOBJS=randomDistTree.o $(DISTTREE_OBJ)
randomDistTree:	$(randomDistTreeOBJS)
	$(CXX) -o $@ $(randomDistTreeOBJS) $(LIBS)
	$(ECHO)

sampleReport.o:  $(COMMON_HPP) 
# $(PHYL_DIR)/distTree.hpp
sampleReportOBJS=sampleReport.o $(CPP_DIR)/common.o
# $(DISTTREE_OBJ)
sampleReport:	$(sampleReportOBJS)
	$(CXX) -o $@ $(sampleReportOBJS) $(LIBS)
	$(ECHO)

statDistTree.o:  $(COMMON_HPP) $(DM_DIR)/dataset.hpp $(PHYL_DIR)/distTree.hpp 
statDistTreeOBJS=statDistTree.o $(DISTTREE_OBJ)
statDistTree:	$(statDistTreeOBJS)
	$(CXX) -o $@ $(statDistTreeOBJS) $(LIBS)
	$(ECHO)

testDistTree.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp 
testDistTreeOBJS=testDistTree.o $(DISTTREE_OBJ)
testDistTree:	$(testDistTreeOBJS)
	$(CXX) -o $@ $(testDistTreeOBJS) $(LIBS)
	$(ECHO)

tree2clusters.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp 
tree2clustersOBJS=tree2clusters.o $(DISTTREE_OBJ)
tree2clusters:	$(tree2clustersOBJS)
	$(CXX) -o $@ $(tree2clustersOBJS) $(LIBS)
	$(ECHO)

tree2genogroup.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp 
tree2genogroupOBJS=tree2genogroup.o $(DISTTREE_OBJ)
tree2genogroup:	$(tree2genogroupOBJS)
	$(CXX) -o $@ $(tree2genogroupOBJS) $(LIBS)
	$(ECHO)

