CODE_DIR=$(shell pwd)/../..
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
  conv_comb \
	dna_diff \
	dna_pair2dissim \
	fasta2dissim \
  feature2dissim \
	feature_request2dissim \
  hash2dissim \
  hash_request2dissim \
  loci_request2dissim \
	main_ortholog \
	mlst2dissim \
	prot_collection2dissim \
	prot_mat2dissim \
	seq2dissim \
	symbet \
	symbet_blastp


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

conv_comb.o:  $(NUMERIC_HPP) 
conv_combOBJS=conv_comb.o $(CPP_DIR)/common.o $(DM_DIR)/numeric.o
conv_comb:	$(conv_combOBJS)
	$(CXX) -o $@ $(conv_combOBJS) $(LIBS)
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

main_ortholog.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/align.hpp
main_orthologOBJS=main_ortholog.o $(ALIGN_OBJ) 
main_ortholog:	$(main_orthologOBJS)
	$(CXX) -o $@ $(main_orthologOBJS) $(LIBS) 
	$(ECHO)

mlst2dissim.o:  $(DM_HPP) 
mlst2dissimOBJS=mlst2dissim.o $(DM_OBJ) 
mlst2dissim:	$(mlst2dissimOBJS)
	$(CXX) -o $@ $(mlst2dissimOBJS) $(LIBS)
	$(ECHO)

prot_collection2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/align.hpp $(DM_DIR)/dataset.hpp $(DISSIM_DIR)/evolution.hpp
prot_collection2dissimOBJS=prot_collection2dissim.o $(ALIGN_OBJ) $(DM_DIR)/dataset.o 
prot_collection2dissim:	$(prot_collection2dissimOBJS)
	$(CXX) -o $@ $(prot_collection2dissimOBJS) $(LIBS) 
	$(ECHO)

prot_mat2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot_mat2dissimOBJS=prot_mat2dissim.o $(SEQ_OBJ)
prot_mat2dissim:	$(prot_mat2dissimOBJS)
	$(CXX) -o $@ $(prot_mat2dissimOBJS) $(LIBS) 
	$(ECHO)

seq2dissim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/align.hpp 
seq2dissimOBJS=seq2dissim.o $(ALIGN_OBJ)
seq2dissim:	$(seq2dissimOBJS)
	$(CXX) -o $@ $(seq2dissimOBJS) $(LIBS) 
	$(ECHO)

symbet.o:  $(NUMERIC_HPP) $(GEN_DIR)/seq.hpp $(DISSIM_DIR)/evolution.hpp
symbetOBJS=symbet.o $(EVOLUTION_OBJ) $(GEN_DIR)/seq.o
symbet:	$(symbetOBJS)
	$(CXX) -o $@ $(symbetOBJS) $(LIBS) 
	$(ECHO)

symbet_blastp.o:  $(NUMERIC_HPP) $(DISSIM_DIR)/evolution.hpp
symbet_blastpOBJS=symbet_blastp.o $(EVOLUTION_OBJ)
symbet_blastp:	$(symbet_blastpOBJS)
	$(CXX) -o $@ $(symbet_blastpOBJS) $(LIBS)
	$(ECHO)

