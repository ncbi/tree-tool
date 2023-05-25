include MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif




############################### Programs ################################

all:	\
  colors \
  curl_easy_test \
  connectPairs \
  effectiveSize \
  extractPairs \
  file2hash \
  index_find \
  list2pairs \
  mergePairs \
  random_words \
  replace_dict \
  setMinus \
  setRandOrd \
  splitList \
  str2hash \
  trav \
  unCgi



colors.o: $(COMMON_HPP)  
colorsOBJS=colors.o $(CPP_DIR)/common.o
colors: $(colorsOBJS)
	$(CXX) -o $@ $(colorsOBJS) $(LIBS)
	$(ECHO)

connectPairs.o: $(COMMON_HPP)  
connectPairsOBJS=connectPairs.o $(CPP_DIR)/common.o
connectPairs: $(connectPairsOBJS)
	$(CXX) -o $@ $(connectPairsOBJS) $(LIBS)
	$(ECHO)

cpp_test.o: $(COMMON_HPP)  
cpp_testOBJS=cpp_test.o $(CPP_DIR)/common.o $(CPP_DIR)/ncurses.o
cpp_test: $(cpp_testOBJS)
	$(CXX) -o $@ $(cpp_testOBJS) $(LIBS) -lncurses  
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

unCgi.o:  $(COMMON_HPP)  
unCgiOBJS=unCgi.o $(CPP_DIR)/common.o
unCgi:	$(unCgiOBJS)
	$(CXX) -o $@ $(unCgiOBJS) $(LIBS)
	$(ECHO)

