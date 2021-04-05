include MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif





############################### Programs ################################

all:	\
  csv2tab \
  effectiveSize \
  extractPairs \
  file2hash \
  group \
  list2pairs \
  mergePairs \
  selColumn \
  setMinus \
  setRandOrd \
  splitList \
  str2hash \
  textTab \
  trav \
  tsv_comp \
  unCgi
	

cpp_test.o:  $(COMMON_HPP)  
cpp_testOBJS=cpp_test.o $(CPP_DIR)/common.o
cpp_test:	$(cpp_testOBJS)
	$(CXX) -o $@ $(cpp_testOBJS) $(LIBS)
	$(ECHO)

csv2tab.o:  $(COMMON_HPP)  
csv2tabOBJS=csv2tab.o $(CPP_DIR)/common.o
csv2tab:	$(csv2tabOBJS)
	$(CXX) -o $@ $(csv2tabOBJS) $(LIBS)
	$(ECHO)

effectiveSize.o:  $(COMMON_HPP)  
effectiveSizeOBJS=effectiveSize.o $(CPP_DIR)/common.o
effectiveSize:	$(effectiveSizeOBJS)
	$(CXX) -o $@ $(effectiveSizeOBJS) $(LIBS)
	$(ECHO)

extractPairs.o:  $(COMMON_HPP)  
extractPairsOBJS=extractPairs.o $(CPP_DIR)/common.o
extractPairs:	$(extractPairsOBJS)
	$(CXX) -o $@ $(extractPairsOBJS) $(LIBS)
	$(ECHO)

file2hash.o:  $(COMMON_HPP)  
file2hashOBJS=file2hash.o $(CPP_DIR)/common.o
file2hash:	$(file2hashOBJS)
	$(CXX) -o $@ $(file2hashOBJS) $(LIBS)
	$(ECHO)

graph_test.o:  $(COMMON_HPP) $(CPP_DIR)/graph.hpp
graph_testOBJS=graph_test.o $(CPP_DIR)/common.o $(CPP_DIR)/graph.o
graph_test:	$(graph_testOBJS)
	$(CXX) -o $@ $(graph_testOBJS) $(LIBS)
	$(ECHO)

group.o:  $(COMMON_HPP)  
groupOBJS=group.o $(CPP_DIR)/common.o
group:	$(groupOBJS)
	$(CXX) -o $@ $(groupOBJS) $(LIBS)
	$(ECHO)

list2pairs.o:  $(COMMON_HPP)  
list2pairsOBJS=list2pairs.o $(CPP_DIR)/common.o
list2pairs:	$(list2pairsOBJS)
	$(CXX) -o $@ $(list2pairsOBJS) $(LIBS)
	$(ECHO)

mergePairs.o:  $(COMMON_HPP)  
mergePairsOBJS=mergePairs.o $(CPP_DIR)/common.o
mergePairs:	$(mergePairsOBJS)
	$(CXX) -o $@ $(mergePairsOBJS) $(LIBS)
	$(ECHO)

selColumn.o:  $(COMMON_HPP)  
selColumnOBJS=selColumn.o $(CPP_DIR)/common.o
selColumn:	$(selColumnOBJS)
	$(CXX) -o $@ $(selColumnOBJS) $(LIBS)
	$(ECHO)

setMinus.o:  $(COMMON_HPP)  
setMinusOBJS=setMinus.o $(CPP_DIR)/common.o
setMinus:	$(setMinusOBJS)
	$(CXX) -o $@ $(setMinusOBJS) $(LIBS)
	$(ECHO)

setRandOrd.o:  $(COMMON_HPP)  
setRandOrdOBJS=setRandOrd.o $(CPP_DIR)/common.o
setRandOrd:	$(setRandOrdOBJS)
	$(CXX) -o $@ $(setRandOrdOBJS) $(LIBS)
	$(ECHO)

splitList.o:  $(COMMON_HPP)  
splitListOBJS=splitList.o $(CPP_DIR)/common.o
splitList:	$(splitListOBJS)
	$(CXX) -o $@ $(splitListOBJS) $(LIBS)
	$(ECHO)

str2hash.o:  $(COMMON_HPP)  
str2hashOBJS=str2hash.o $(CPP_DIR)/common.o
str2hash:	$(str2hashOBJS)
	$(CXX) -o $@ $(str2hashOBJS) $(LIBS)
	$(ECHO)

textTab.o:  $(COMMON_HPP)  
textTabOBJS=textTab.o $(CPP_DIR)/common.o
textTab:	$(textTabOBJS)
	$(CXX) -o $@ $(textTabOBJS) $(LIBS)
	$(ECHO)

trav.o:  $(COMMON_HPP)  
travOBJS=trav.o $(CPP_DIR)/common.o
trav:	$(travOBJS)
	$(CXX) -o $@ $(travOBJS) $(LIBS)
	$(ECHO)

tsv_comp.o:  $(COMMON_HPP)  
tsv_compOBJS=tsv_comp.o $(CPP_DIR)/common.o
tsv_comp:	$(tsv_compOBJS)
	$(CXX) -o $@ $(tsv_compOBJS) $(LIBS)
	$(ECHO)

unCgi.o:  $(COMMON_HPP)  
unCgiOBJS=unCgi.o $(CPP_DIR)/common.o
unCgi:	$(unCgiOBJS)
	$(CXX) -o $@ $(unCgiOBJS) $(LIBS)
	$(ECHO)

