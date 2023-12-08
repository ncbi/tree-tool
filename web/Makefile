include ../MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif




############################### Programs ################################

all:  unCgi


unCgi.o:  $(COMMON_HPP)  
unCgiOBJS=unCgi.o $(WEB_OBJ)
unCgi:	$(unCgiOBJS)
	$(CXX) -o $@ $(unCgiOBJS) $(LIBS)
	$(ECHO)

