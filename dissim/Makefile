include ../MakeRules
CPPFLAGS=$(OptCPPFLAGS)
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
  additive \
	blast2ani \
	blast2cons \
  combine_dissims \
	dna2dissim \
	dna_diff \
	dna_pair2dissim \
	fasta2dissim \
  feature2dissim \
	feature_request2dissim \
  hash2dissim \
  hash_request2dissim \
  loci_request2dissim \
	mlst2dissim \
	prot2dissim \
	prots_pair2dissim


additive.o:  $(DM_HPP) 
additiveOBJS=additive.o $(DM_OBJ) 
additive:	$(additiveOBJS)
	$(CXX) -o $@ $(additiveOBJS) $(LIBS)
	$(ECHO)

blast2ani.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp
blast2aniOBJS=blast2ani.o $(NUMERIC_OBJ)
blast2ani:	$(blast2aniOBJS)
	$(CXX) -o $@ $(blast2aniOBJS) $(LIBS) 
	$(ECHO)

blast2consOBJS=blast2cons.o
blast2cons:	$(blast2consOBJS)
	$(CXX) -o $@ $(blast2consOBJS) $(LIBS) 
	$(ECHO)

combine_dissims.o:  $(NUMERIC_HPP) 
combine_dissimsOBJS=combine_dissims.o $(CPP_DIR)/common.o $(DM_DIR)/numeric.o
combine_dissims:	$(combine_dissimsOBJS)
	$(CXX) -o $@ $(combine_dissimsOBJS) $(LIBS)
	$(ECHO)

dna2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/align.hpp 
dna2dissimOBJS=dna2dissim.o $(ALIGN_OBJ)
dna2dissim:	$(dna2dissimOBJS)
	$(CXX) -o $@ $(dna2dissimOBJS) $(LIBS) 
	$(ECHO)

dna_diff.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_diffOBJS=dna_diff.o $(SEQ_OBJ)
dna_diff:	$(dna_diffOBJS)
	$(CXX) -o $@ $(dna_diffOBJS) $(LIBS) 
	$(ECHO)

dna_pair2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/align.hpp 
dna_pair2dissimOBJS=dna_pair2dissim.o $(ALIGN_OBJ) 
dna_pair2dissim:	$(dna_pair2dissimOBJS)
	$(CXX) -o $@ $(dna_pair2dissimOBJS) $(LIBS) 
	$(ECHO)

fasta2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/align.hpp $(DM_DIR)/dataset.hpp $(DM_DIR)/matrix.hpp 
fasta2dissimOBJS=fasta2dissim.o $(ALIGN_OBJ) $(DM_DIR)/dataset.o 
fasta2dissim:	$(fasta2dissimOBJS)
	$(CXX) -o $@ $(fasta2dissimOBJS) $(LIBS) 
	$(ECHO)

feature2dissim.o:  $(DM_HPP) $(DISSIM_DIR)/evolution.hpp
feature2dissimOBJS=feature2dissim.o $(DM_OBJ) $(DISSIM_DIR)/evolution.o
feature2dissim:	$(feature2dissimOBJS)
	$(CXX) -o $@ $(feature2dissimOBJS) $(LIBS)
	$(ECHO)

feature_request2dissim.o:  $(NUMERIC_HPP) $(DISSIM_DIR)/evolution.hpp
feature_request2dissimOBJS=feature_request2dissim.o $(DM_OBJ) $(DISSIM_DIR)/evolution.o
feature_request2dissim:	$(feature_request2dissimOBJS)
	$(CXX) -o $@ $(feature_request2dissimOBJS) $(LIBS)
	$(ECHO)

hash2dissim.o:  $(DM_HPP) $(DISSIM_DIR)/evolution.hpp 
hash2dissimOBJS=hash2dissim.o $(DM_OBJ) $(DISSIM_DIR)/evolution.o
hash2dissim:	$(hash2dissimOBJS)
	$(CXX) -o $@ $(hash2dissimOBJS) $(LIBS)
	$(ECHO)

hash_request2dissim.o:  $(NUMERIC_HPP) $(DISSIM_DIR)/evolution.hpp 
hash_request2dissimOBJS=hash_request2dissim.o $(DM_OBJ) $(DISSIM_DIR)/evolution.o
hash_request2dissim:	$(hash_request2dissimOBJS)
	$(CXX) -o $@ $(hash_request2dissimOBJS) $(LIBS)
	$(ECHO)

loci_request2dissim.o:  $(COMMON_HPP) 
loci_request2dissimOBJS=loci_request2dissim.o $(CPP_DIR)/common.o
loci_request2dissim:	$(loci_request2dissimOBJS)
	$(CXX) -o $@ $(loci_request2dissimOBJS) $(LIBS)
	$(ECHO)

mlst2dissim.o:  $(DM_HPP) 
mlst2dissimOBJS=mlst2dissim.o $(DM_OBJ) 
mlst2dissim:	$(mlst2dissimOBJS)
	$(CXX) -o $@ $(mlst2dissimOBJS) $(LIBS)
	$(ECHO)

prot2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot2dissimOBJS=prot2dissim.o $(SEQ_OBJ)
prot2dissim:	$(prot2dissimOBJS)
	$(CXX) -o $@ $(prot2dissimOBJS) $(LIBS) 
	$(ECHO)

prots_pair2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/align.hpp $(DM_DIR)/dataset.hpp $(DISSIM_DIR)/evolution.hpp
prots_pair2dissimOBJS=prots_pair2dissim.o $(ALIGN_OBJ) $(DM_DIR)/dataset.o 
prots_pair2dissim:	$(prots_pair2dissimOBJS)
	$(CXX) -o $@ $(prots_pair2dissimOBJS) $(LIBS) 
	$(ECHO)

