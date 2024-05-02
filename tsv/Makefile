CODE_DIR=$(shell pwd)/../..
include ../MakeRules
CPPFLAGS = $(OptCPPFLAGS) 
#CPPFLAGS = $(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif




############################### Programs ################################

ALL=	\
  csv2tab \
  fixed2tsv \
  pairs2tsv \
  triple2tsv \
  tsv2insert \
  tsv2triple \
  tsv_aggr_comp \
  tsv_cat \
  tsv_cluster \
  tsv_comp \
  tsv_group \
  tsv_join \
  tsv_rename \
  tsv_schema \
  tsv_split

ifdef AT_NCBI
  all: $(ALL) \
    tsv_view
else
  all:  $(ALL)
endif



csv2tab.o:  $(COMMON_HPP)  
csv2tabOBJS=csv2tab.o $(CPP_DIR)/common.o
csv2tab:  $(csv2tabOBJS)
	$(CXX) -o $@ $(csv2tabOBJS) $(LIBS)
	$(ECHO)

fixed2tsv.o:  $(COMMON_HPP) 
fixed2tsvOBJS=fixed2tsv.o $(CPP_DIR)/common.o 
fixed2tsv:	$(fixed2tsvOBJS)
	$(CXX) -o $@ $(fixed2tsvOBJS) $(LIBS)
	$(ECHO)

pairs2tsv.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
pairs2tsvOBJS=pairs2tsv.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
pairs2tsv:	$(pairs2tsvOBJS)
	$(CXX) -o $@ $(pairs2tsvOBJS) $(LIBS)
	$(ECHO)

triple2tsv.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
triple2tsvOBJS=triple2tsv.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
triple2tsv:	$(triple2tsvOBJS)
	$(CXX) -o $@ $(triple2tsvOBJS) $(LIBS)
	$(ECHO)

tsv2insert.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv2insertOBJS=tsv2insert.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv2insert:	$(tsv2insertOBJS)
	$(CXX) -o $@ $(tsv2insertOBJS) $(LIBS)
	$(ECHO)

tsv2triple.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv2tripleOBJS=tsv2triple.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv2triple:	$(tsv2tripleOBJS)
	$(CXX) -o $@ $(tsv2tripleOBJS) $(LIBS)
	$(ECHO)

tsv_aggr_comp.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv_aggr_compOBJS=tsv_aggr_comp.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv_aggr_comp:	$(tsv_aggr_compOBJS)
	$(CXX) -o $@ $(tsv_aggr_compOBJS) $(LIBS)
	$(ECHO)

tsv_cat.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv_catOBJS=tsv_cat.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv_cat:	$(tsv_catOBJS)
	$(CXX) -o $@ $(tsv_catOBJS) $(LIBS)
	$(ECHO)

tsv_cluster.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp $(DM_DIR)/numeric.hpp $(DISSIM_DIR)/evolution.hpp
tsv_clusterOBJS=tsv_cluster.o $(TSV_DIR)/tsv.o $(DISSIM_DIR)/evolution.o $(DM_OBJ)
tsv_cluster:	$(tsv_clusterOBJS)
	$(CXX) -o $@ $(tsv_clusterOBJS) $(LIBS)
	$(ECHO)

tsv_comp.o:  $(COMMON_HPP) 
tsv_compOBJS=tsv_comp.o $(CPP_DIR)/common.o 
tsv_comp:	$(tsv_compOBJS)
	$(CXX) -o $@ $(tsv_compOBJS) $(LIBS)
	$(ECHO)

tsv_group.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv_groupOBJS=tsv_group.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv_group:	$(tsv_groupOBJS)
	$(CXX) -o $@ $(tsv_groupOBJS) $(LIBS)
	$(ECHO)

tsv_join.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv_joinOBJS=tsv_join.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv_join: $(tsv_joinOBJS)
	$(CXX) -o $@ $(tsv_joinOBJS) $(LIBS)
	$(ECHO)

tsv_rename.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv_renameOBJS=tsv_rename.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv_rename:	$(tsv_renameOBJS)
	$(CXX) -o $@ $(tsv_renameOBJS) $(LIBS)
	$(ECHO)

tsv_schema.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv_schemaOBJS=tsv_schema.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv_schema: $(tsv_schemaOBJS)
	$(CXX) -o $@ $(tsv_schemaOBJS) $(LIBS)
	$(ECHO)

tsv_split.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
tsv_splitOBJS=tsv_split.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
tsv_split:	$(tsv_splitOBJS)
	$(CXX) -o $@ $(tsv_splitOBJS) $(LIBS)
	$(ECHO)

tsv_view.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp $(CPP_DIR)/ncurses.hpp 
tsv_viewOBJS=tsv_view.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o $(CPP_DIR)/ncurses.o 
tsv_view:	$(tsv_viewOBJS)
	$(CXX) -o $@ $(tsv_viewOBJS) $(LIBS) -lncurses
	$(ECHO)

