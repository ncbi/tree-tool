include $(HOME)/code/cpp/MakeRules
CPPFLAGS=$(OptCPPFLAGS) 
#CPPFLAGS=$(DebugCPPFLAGS)




############################### Programs ################################

all:	\
  grammar_test
	

grammar_test.o:  $(COMMON_HPP) $(CPP_DIR)/graph.hpp cfgrammar.hpp
grammar_testOBJS=grammar_test.o $(CPP_DIR)/common.o $(CPP_DIR)/graph.o $(LANG_DIR)/cfgrammar.o
grammar_test:	$(grammar_testOBJS)
	$(CXX) -o $@ $(grammar_testOBJS) $(LIBS)
	$(ECHO)

