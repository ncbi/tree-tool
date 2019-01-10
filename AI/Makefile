include $(HOME)/code/cpp/MakeRules
CPPFLAGS=$(OptCPPFLAGS) 
#CPPFLAGS=$(DebugCPPFLAGS)




############################### Programs ################################

all:  \
  algebra
	

algebra.o:  $(COMMON_HPP) $(AI_DIR)/expr.hpp
algebraOBJS=algebra.o $(CPP_DIR)/common.o $(AI_DIR)/expr.o
algebra:	$(algebraOBJS)
	$(CXX) -o $@ $(algebraOBJS) $(LIBS)
	$(ECHO)

