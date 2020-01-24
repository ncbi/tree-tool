include ../MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif





############################### Programs ################################

all:	\
  xml2schema \
  xml_merge_schemas \
  xml_schema2ddl \
  xml_schema2flat
	

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

xml_schema2ddl.o:  $(COMMON_HPP) $(XML_DIR)/xml.hpp
xml_schema2ddlOBJS=xml_schema2ddl.o $(CPP_DIR)/common.o $(XML_DIR)/xml.o
xml_schema2ddl:	$(xml_schema2ddlOBJS)
	$(CXX) -o $@ $(xml_schema2ddlOBJS) $(LIBS)
	$(ECHO)

xml_schema2flat.o:  $(COMMON_HPP) $(XML_DIR)/xml.hpp
xml_schema2flatOBJS=xml_schema2flat.o $(CPP_DIR)/common.o $(XML_DIR)/xml.o
xml_schema2flat:	$(xml_schema2flatOBJS)
	$(CXX) -o $@ $(xml_schema2flatOBJS) $(LIBS)
	$(ECHO)

