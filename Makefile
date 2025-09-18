CODE_DIR=$(shell pwd)/..
include MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif




############################### Programs ################################

ALL=	\
  ascii \
  colors_test \
  connectPairs \
  curl_easy_test \
  effectiveSize \
  extractPairs \
  file2hash \
  graph_test \
  gzip_test \
  index_find \
  key_test \
  list2pairs \
  mergePairs \
  min_spanning_forest \
  multilist2subset \
  objHash_find \
  random_words \
  replace_dict \
  setMinus \
  setRandOrd \
  splitList \
  str2hash \
  trav

ifdef AT_NCBI
  all:	$(ALL) \
    ncurses_key_test
else
  all:	$(ALL)
endif


ascii.o: $(COMMON_HPP)  
asciiOBJS=ascii.o $(CPP_DIR)/common.o
ascii: $(asciiOBJS)
	$(CXX) -o $@ $(asciiOBJS) $(LIBS)
	$(ECHO)

colors_test.o: $(COMMON_HPP)  
colors_testOBJS=colors_test.o $(CPP_DIR)/common.o
colors_test: $(colors_testOBJS)
	$(CXX) -o $@ $(colors_testOBJS) $(LIBS)
	$(ECHO)

connectPairs.o: $(COMMON_HPP)  
connectPairsOBJS=connectPairs.o $(CPP_DIR)/common.o
connectPairs: $(connectPairsOBJS)
	$(CXX) -o $@ $(connectPairsOBJS) $(LIBS)
	$(ECHO)

cpp_test.o: $(COMMON_HPP) $(XML_DIR)/xml.hpp 
cpp_testOBJS=cpp_test.o $(CPP_DIR)/common.o $(XML_DIR)/xml.o $(TSV_DIR)/tsv.o
cpp_test: $(cpp_testOBJS)
	$(CXX) -o $@ $(cpp_testOBJS) $(LIBS) 
	  # -lncursesw -lz
	$(ECHO)

curl_easy_test.o: $(COMMON_HPP) $(CPP_DIR)/curl_easy.hpp
curl_easy_testOBJS=curl_easy_test.o $(CPP_DIR)/common.o $(CPP_DIR)/curl_easy.o
curl_easy_test: $(curl_easy_testOBJS)
	$(CXX) -o $@ $(curl_easy_testOBJS) $(LIBS) -lcurl
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

fstream_test.o: $(COMMON_HPP)  
fstream_testOBJS=fstream_test.o $(CPP_DIR)/common.o
fstream_test: $(fstream_testOBJS)
	$(CXX) -o $@ $(fstream_testOBJS) $(LIBS)
	$(ECHO)

graph_test.o:  $(COMMON_HPP) $(CPP_DIR)/graph.hpp
graph_testOBJS=graph_test.o $(CPP_DIR)/common.o $(CPP_DIR)/graph.o
graph_test:	$(graph_testOBJS)
	$(CXX) -o $@ $(graph_testOBJS) $(LIBS)
	$(ECHO)

gzip_test.o: $(COMMON_HPP) $(CPP_DIR)/gzip.hpp
gzip_testOBJS=gzip_test.o $(CPP_DIR)/common.o $(CPP_DIR)/gzip.o 
gzip_test: $(gzip_testOBJS)
	$(CXX) -o $@ $(gzip_testOBJS) $(LIBS) -lz
	$(ECHO)

index_find.o:  $(COMMON_HPP)  
index_findOBJS=index_find.o $(CPP_DIR)/common.o
index_find:	$(index_findOBJS)
	$(CXX) -o $@ $(index_findOBJS) $(LIBS)
	$(ECHO)

key_test.o: $(COMMON_HPP)  
key_testOBJS=key_test.o $(CPP_DIR)/common.o
key_test: $(key_testOBJS)
	$(CXX) -o $@ $(key_testOBJS) $(LIBS)
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

min_spanning_forest.o: $(COMMON_HPP) $(CPP_DIR)/graph.hpp
min_spanning_forestOBJS=min_spanning_forest.o $(CPP_DIR)/common.o $(CPP_DIR)/graph.o
min_spanning_forest: $(min_spanning_forestOBJS)
	$(CXX) -o $@ $(min_spanning_forestOBJS) $(LIBS)
	$(ECHO)

multilist2subset.o:  $(COMMON_HPP)  
multilist2subsetOBJS=multilist2subset.o $(CPP_DIR)/common.o
multilist2subset:  $(multilist2subsetOBJS)
	$(CXX) -o $@ $(multilist2subsetOBJS) $(LIBS)
	$(ECHO)

ncurses_key_test.o: $(COMMON_HPP) $(CPP_DIR)/ncurses.hpp
ncurses_key_testOBJS=ncurses_key_test.o $(CPP_DIR)/common.o $(CPP_DIR)/ncurses.o
ncurses_key_test: $(ncurses_key_testOBJS)
	$(CXX) -o $@ $(ncurses_key_testOBJS) $(LIBS) -lncursesw
	$(ECHO)

objHash_find.o:  $(COMMON_HPP) 
objHash_findOBJS=objHash_find.o $(CPP_DIR)/common.o
objHash_find:	$(objHash_findOBJS)
	$(CXX) -o $@ $(objHash_findOBJS) $(LIBS) 
	$(ECHO)

random_words.o: $(COMMON_HPP)  
random_wordsOBJS=random_words.o $(CPP_DIR)/common.o
random_words: $(random_wordsOBJS)
	$(CXX) -o $@ $(random_wordsOBJS) $(LIBS)
	$(ECHO)

replace_dict.o: $(COMMON_HPP)  
replace_dictOBJS=replace_dict.o $(CPP_DIR)/common.o
replace_dict: $(replace_dictOBJS)
	$(CXX) -o $@ $(replace_dictOBJS) $(LIBS)
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

