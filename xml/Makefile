include ../MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif





############################### Programs ################################

all:	\
  xml2schema \
  xml_merge_schemas
	

xml2schema.o:  $(COMMON_HPP) $(XML_DIR)/xml.hpp
xml2schemaOBJS=xml2schema.o $(CPP_DIR)/common.o $(XML_DIR)/xml.o
xml2schema:	$(xml2schemaOBJS)
	$(CXX) -o $@ $(xml2schemaOBJS) $(LIBS)
	$(ECHO)

xml_merge_schemas.o:  $(COMMON_HPP) $(XML_DIR)/xml.hpp
xml_merge_schemasOBJS=xml_merge_schemas.o $(CPP_DIR)/common.o $(XML_DIR)/xml.o
xml_merge_schemas:	$(xml_merge_schemasOBJS)
	$(CXX) -o $@ $(xml_merge_schemasOBJS) $(LIBS)
	$(ECHO)

