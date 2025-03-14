CODE_DIR=$(shell pwd)/../..
include ../MakeRules
CPPFLAGS=$(OptCPPFLAGS) 
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
	asnt2tree \
  compareTrees \
  distTree_new \
  distTree_refresh_dissim \
  dm2feature \
	fasta2feature \
  feature2gain_loss \
	makeDistTree \
	makeDistTree-dev \
	makeFeatureTree \
	mutation2feature \
	newick2tree \
	printDistTree \
	randomDistTree \
	replaceDistTree_match \
	replaceDistTree_reroot \
	statDistTree \
	testDistTree \
	tree2genogroup \
	tree2indiscern
	

asnt2tree.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp $(PHYL_DIR)/distTree.hpp 
asnt2treeOBJS=asnt2tree.o $(DISTTREE_OBJ) 
asnt2tree:	$(asnt2treeOBJS)
	$(CXX) -o $@ $(asnt2treeOBJS) $(LIBS) 
	$(ECHO)

dm2feature.o:  $(DM_HPP) 
dm2featureOBJS=dm2feature.o $(DM_OBJ) 
dm2feature:	$(dm2featureOBJS)
	$(CXX) -o $@ $(dm2featureOBJS) $(LIBS)
	$(ECHO)

compareTrees.o:  $(COMMON_HPP) $(PHYL_DIR)/distTree.hpp $(PHYL_DIR)/featureTree.hpp 
compareTreesOBJS=compareTrees.o $(DISTTREE_OBJ) $(PHYL_DIR)/featureTree.o
compareTrees:	$(compareTreesOBJS)
	$(CXX) -o $@ $(compareTreesOBJS) $(LIBS)
	$(ECHO)
	
distTree_new.o:  $(DISTTREE_HPP)
distTree_newOBJS=distTree_new.o $(DISTTREE_OBJ)
distTree_new:	$(distTree_newOBJS)
	$(CXX) -o $@ $(distTree_newOBJS) $(LIBS)
	$(ECHO)

distTree_refresh_dissim.o:  $(DISTTREE_HPP)
distTree_refresh_dissimOBJS=distTree_refresh_dissim.o $(DISTTREE_OBJ)
distTree_refresh_dissim:	$(distTree_refresh_dissimOBJS)
	$(CXX) -o $@ $(distTree_refresh_dissimOBJS) $(LIBS)
	$(ECHO)

fasta2feature.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
fasta2featureOBJS=fasta2feature.o $(SEQ_OBJ)
fasta2feature:	$(fasta2featureOBJS)
	$(CXX) -o $@ $(fasta2featureOBJS) $(LIBS) 
	$(ECHO)

feature2gain_loss.o:  $(COMMON_HPP) $(PHYL_DIR)/featureTree.hpp 
feature2gain_lossOBJS=feature2gain_loss.o $(FEATURE_TREE_OBJ)
feature2gain_loss:	$(feature2gain_lossOBJS)
	$(CXX) -o $@ $(feature2gain_lossOBJS) $(LIBS)
	$(ECHO)

makeDistTree.o:  $(DISTTREE_HPP)
  # $(COMMON_HPP) $(CPP_DIR)/graph.hpp $(DM_DIR)/dataset.hpp $(PHYL_DIR)/distTree.hpp 
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

mutation2feature.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp
mutation2featureOBJS=mutation2feature.o $(SEQ_OBJ)
mutation2feature:	$(mutation2featureOBJS)
	$(CXX) -o $@ $(mutation2featureOBJS) $(LIBS)
	$(ECHO)

newick2tree.o:  $(DISTTREE_HPP)
newick2treeOBJS=newick2tree.o $(DISTTREE_OBJ)
newick2tree:	$(newick2treeOBJS)
	$(CXX) -o $@ $(newick2treeOBJS) $(LIBS)
	$(ECHO)

printDistTree.o:  $(DISTTREE_HPP)
printDistTreeOBJS=printDistTree.o $(DISTTREE_OBJ)
printDistTree:	$(printDistTreeOBJS)
	$(CXX) -o $@ $(printDistTreeOBJS) $(LIBS)
	$(ECHO)

randomDistTree.o:  $(DISTTREE_HPP)
randomDistTreeOBJS=randomDistTree.o $(DISTTREE_OBJ)
randomDistTree:	$(randomDistTreeOBJS)
	$(CXX) -o $@ $(randomDistTreeOBJS) $(LIBS)
	$(ECHO)

replaceDistTree_match.o:  $(DISTTREE_HPP)
replaceDistTree_matchOBJS=replaceDistTree_match.o $(DISTTREE_OBJ)
replaceDistTree_match:	$(replaceDistTree_matchOBJS)
	$(CXX) -o $@ $(replaceDistTree_matchOBJS) $(LIBS)
	$(ECHO)

replaceDistTree_reroot.o:  $(DISTTREE_HPP)
replaceDistTree_rerootOBJS=replaceDistTree_reroot.o $(DISTTREE_OBJ)
replaceDistTree_reroot:	$(replaceDistTree_rerootOBJS)
	$(CXX) -o $@ $(replaceDistTree_rerootOBJS) $(LIBS)
	$(ECHO)

statDistTree.o:  $(DISTTREE_HPP)
statDistTreeOBJS=statDistTree.o $(DISTTREE_OBJ)
statDistTree:	$(statDistTreeOBJS)
	$(CXX) -o $@ $(statDistTreeOBJS) $(LIBS)
	$(ECHO)

testDistTree.o:  $(DISTTREE_HPP)
testDistTreeOBJS=testDistTree.o $(DISTTREE_OBJ)
testDistTree:	$(testDistTreeOBJS)
	$(CXX) -o $@ $(testDistTreeOBJS) $(LIBS)
	$(ECHO)

tree2genogroup.o:  $(DISTTREE_HPP)
tree2genogroupOBJS=tree2genogroup.o $(DISTTREE_OBJ)
tree2genogroup:	$(tree2genogroupOBJS)
	$(CXX) -o $@ $(tree2genogroupOBJS) $(LIBS)
	$(ECHO)

tree2indiscern.o:  $(DISTTREE_HPP)
tree2indiscernOBJS=tree2indiscern.o $(DISTTREE_OBJ)
tree2indiscern:	$(tree2indiscernOBJS)
	$(CXX) -o $@ $(tree2indiscernOBJS) $(LIBS)
	$(ECHO)

