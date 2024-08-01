CODE_DIR=$(shell pwd)/../../..
include ../../MakeRules
CPPFLAGS=$(OptCPPFLAGS) 
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
  attr2_2pairs \
  attr2_2paup \
  attr2_2phylip \
	obj_attr2dm \
	pairs2dm \
	tsv2dm
	

attr2_2pairs.o:  $(DM_HPP) 
attr2_2pairsOBJS=attr2_2pairs.o $(DM_OBJ) 
attr2_2pairs:	$(attr2_2pairsOBJS)
	$(CXX) -o $@ $(attr2_2pairsOBJS) $(LIBS)
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

obj_attr2dm.o:  $(DM_HPP) 
obj_attr2dmOBJS=obj_attr2dm.o $(DM_OBJ) 
obj_attr2dm:	$(obj_attr2dmOBJS)
	$(CXX) -o $@ $(obj_attr2dmOBJS) $(LIBS)
	$(ECHO)

pairs2dm.o:  $(DM_HPP) 
pairs2dmOBJS=pairs2dm.o $(DM_OBJ) 
pairs2dm:	$(pairs2dmOBJS)
	$(CXX) -o $@ $(pairs2dmOBJS) $(LIBS)
	$(ECHO)

tsv2dm.o:  $(DM_HPP) 
tsv2dmOBJS=tsv2dm.o $(DM_OBJ) $(TSV_DIR)/tsv.o
tsv2dm:	$(tsv2dmOBJS)
	$(CXX) -o $@ $(tsv2dmOBJS) $(LIBS)
	$(ECHO)

