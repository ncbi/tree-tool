include ../MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif




############################### Programs ################################

ALL=	\
  tsv_cat \
  tsv_comp \
  tsv_group \
  tsv_join \
  tsv_schema

ifdef AT_NCBI
  all: $(ALL) \
    tsv_view
else
  all:  $(ALL)
endif



tsv_cat.o:  $(COMMON_HPP)  
tsv_catOBJS=tsv_cat.o $(CPP_DIR)/common.o
tsv_cat:	$(tsv_catOBJS)
	$(CXX) -o $@ $(tsv_catOBJS) $(LIBS)
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

