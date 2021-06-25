include MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif




############################### Programs ################################

ALL=	\
  connectPairs \
  csv2tab \
  effectiveSize \
  extractPairs \
  file2hash \
  index_find \
  list2pairs \
  mergePairs \
  selColumn \
  setMinus \
  setRandOrd \
  splitList \
  str2hash \
  trav \
  tsv_comp \
  tsv_group \
  tsv_join \
  tsv_schema \
  unCgi

ifdef AT_NCBI
  all: $(ALL) \
    tsv_view
else
  all:  $(ALL)
endif



connectPairs.o: $(COMMON_HPP)  
connectPairsOBJS=connectPairs.o $(CPP_DIR)/common.o
connectPairs: $(connectPairsOBJS)
	$(CXX) -o $@ $(connectPairsOBJS) $(LIBS)
	$(ECHO)

cpp_test.o: $(COMMON_HPP)  
cpp_testOBJS=cpp_test.o $(CPP_DIR)/common.o
cpp_test: $(cpp_testOBJS)
	$(CXX) -o $@ $(cpp_testOBJS) $(LIBS)
	$(ECHO)

csv2tab.o:  $(COMMON_HPP)  
csv2tabOBJS=csv2tab.o $(CPP_DIR)/common.o
csv2tab:  $(csv2tabOBJS)
	$(CXX) -o $@ $(csv2tabOBJS) $(LIBS)
	$(ECHO)

effectiveSize.o:  $(COMMON_HPP)  
effectiveSizeOBJS=effectiveSize.o $(CPP_DIR)/common.o
effectiveSize:  $(effectiveSizeOBJS)
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

index_find.o:  $(COMMON_HPP)  
index_findOBJS=index_find.o $(CPP_DIR)/common.o
index_find:	$(index_findOBJS)
	$(CXX) -o $@ $(index_findOBJS) $(LIBS)
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

tsv_group.o:  $(COMMON_HPP)  
tsv_groupOBJS=tsv_group.o $(CPP_DIR)/common.o
tsv_group:	$(tsv_groupOBJS)
	$(CXX) -o $@ $(tsv_groupOBJS) $(LIBS)
	$(ECHO)

tsv_join.o:  $(COMMON_HPP)  
tsv_joinOBJS=tsv_join.o $(CPP_DIR)/common.o
tsv_join: $(tsv_joinOBJS)
	$(CXX) -o $@ $(tsv_joinOBJS) $(LIBS)
	$(ECHO)

tsv_schema.o:  $(COMMON_HPP)  
tsv_schemaOBJS=tsv_schema.o $(CPP_DIR)/common.o
tsv_schema: $(tsv_schemaOBJS)
	$(CXX) -o $@ $(tsv_schemaOBJS) $(LIBS)
	$(ECHO)

tsv_view.o:  $(COMMON_HPP) $(CPP_DIR)/ncurses.hpp 
tsv_viewOBJS=tsv_view.o $(CPP_DIR)/common.o $(CPP_DIR)/ncurses.o 
tsv_view:	$(tsv_viewOBJS)
	$(CXX) -o $@ $(tsv_viewOBJS) $(LIBS) -lncurses
	$(ECHO)

unCgi.o:  $(COMMON_HPP)  
unCgiOBJS=unCgi.o $(CPP_DIR)/common.o
unCgi:	$(unCgiOBJS)
	$(CXX) -o $@ $(unCgiOBJS) $(LIBS)
	$(ECHO)

