include ../MakeRules
CPPFLAGS=$(OptCPPFLAGS)
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
	asn_grep \
	blastmat \
	blastp_merge \
	dna2orfs \
	dna2stat \
	dna_coverage \
	dna_cut \
	dna_rev \
	dna_trim \
	extractFastaDna \
	extractFastaProt \
	fasta2GC \
	fasta2hash \
	fasta2len \
	fastaAddGi \
	fastaDna2pairs \
	findProt \
	GeneMark2CDS \
	hmm2prot \
	hmm_tc1 \
	hmmAddCutoff \
	hmmNAME2ACC \
	hmmsearch2besthits \
	kmerIndex_add \
	kmerIndex_find \
	kmerIndex_make \
	kmerIndex_stat \
	mlst2hash \
	mutation_tab \
	orf2prot \
	prot2fingerprints \
	prot2triplets \
	prot_check \
	prot_fingerprint_grep \
	prot_grep_short \
	prots2hmm_signature \
	prots_pair2stat \
	replaceFastaHeader \
	replaceFastaNames \
	splitFasta \
	splitHmmLib \
	tblastn2orfs \
	uniqSeq \
	uniqProtRef


asn_grep.o:	$(COMMON_HPP) $(GEN_DIR)/genbank_text.hpp
asn_grepOBJS=asn_grep.o $(CPP_DIR)/common.o $(GEN_DIR)/genbank_text.o
asn_grep:	$(asn_grepOBJS)
	$(CXX) -o $@ $(asn_grepOBJS) $(LIBS) 
	$(ECHO)

blastmat.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
blastmatOBJS=blastmat.o $(SEQ_OBJ)
blastmat:	$(blastmatOBJS)
	$(CXX) -o $@ $(blastmatOBJS) $(LIBS) 
	$(ECHO)

blastp_merge.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp $(DM_DIR)/matrix.hpp
blastp_mergeOBJS=blastp_merge.o $(CPP_DIR)/common.o $(DM_DIR)/numeric.o $(DM_DIR)/matrix.o
blastp_merge:	$(blastp_mergeOBJS)
	$(CXX) -o $@ $(blastp_mergeOBJS) $(LIBS) 
	$(ECHO)

dna2orfs.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna2orfsOBJS=dna2orfs.o $(SEQ_OBJ)
dna2orfs:	$(dna2orfsOBJS)
	$(CXX) -o $@ $(dna2orfsOBJS) $(LIBS) 
	$(ECHO)

dna2stat.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna2statOBJS=dna2stat.o $(SEQ_OBJ)
dna2stat:	$(dna2statOBJS)
	$(CXX) -o $@ $(dna2statOBJS) $(LIBS) 
	$(ECHO)

dna_coverage.o:  $(COMMON_HPP) 
dna_coverageOBJS=dna_coverage.o $(CPP_DIR)/common.o
dna_coverage:	$(dna_coverageOBJS)
	$(CXX) -o $@ $(dna_coverageOBJS) $(LIBS) 
	$(ECHO)

dna_cut.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_cutOBJS=dna_cut.o $(SEQ_OBJ)
dna_cut:	$(dna_cutOBJS)
	$(CXX) -o $@ $(dna_cutOBJS) $(LIBS) 
	$(ECHO)

dna_rev.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_revOBJS=dna_rev.o $(SEQ_OBJ)
dna_rev:	$(dna_revOBJS)
	$(CXX) -o $@ $(dna_revOBJS) $(LIBS) 
	$(ECHO)

dna_trim.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_trimOBJS=dna_trim.o $(SEQ_OBJ)
dna_trim:	$(dna_trimOBJS)
	$(CXX) -o $@ $(dna_trimOBJS) $(LIBS) 
	$(ECHO)

extractFastaDna.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
extractFastaDnaOBJS=extractFastaDna.o $(SEQ_OBJ)
extractFastaDna:	$(extractFastaDnaOBJS)
	$(CXX) -o $@ $(extractFastaDnaOBJS) $(LIBS) 
	$(ECHO)

extractFastaProt.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
extractFastaProtOBJS=extractFastaProt.o $(SEQ_OBJ)
extractFastaProt:	$(extractFastaProtOBJS)
	$(CXX) -o $@ $(extractFastaProtOBJS) $(LIBS) 
	$(ECHO)

fasta2GC.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DM_HPP) 
fasta2GCOBJS=fasta2GC.o $(GEN_DIR)/seq.o $(DM_OBJ)
fasta2GC:	$(fasta2GCOBJS)
	$(CXX) -o $@ $(fasta2GCOBJS) $(LIBS) 
	$(ECHO)

fasta2hash.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
fasta2hashOBJS=fasta2hash.o $(SEQ_OBJ)
fasta2hash:	$(fasta2hashOBJS)
	$(CXX) -o $@ $(fasta2hashOBJS) $(LIBS) 
	$(ECHO)

fasta2len.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
fasta2lenOBJS=fasta2len.o $(SEQ_OBJ)
fasta2len:	$(fasta2lenOBJS)
	$(CXX) -o $@ $(fasta2lenOBJS) $(LIBS) 
	$(ECHO)

fastaAddGi.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
fastaAddGiOBJS=fastaAddGi.o $(SEQ_OBJ)
fastaAddGi:	$(fastaAddGiOBJS)
	$(CXX) -o $@ $(fastaAddGiOBJS) $(LIBS) 
	$(ECHO)

fastaDna2pairs.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
fastaDna2pairsOBJS=fastaDna2pairs.o $(SEQ_OBJ)
fastaDna2pairs:	$(fastaDna2pairsOBJS)
	$(CXX) -o $@ $(fastaDna2pairsOBJS) $(LIBS) 
	$(ECHO)

findProt.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
findProtOBJS=findProt.o $(SEQ_OBJ)
findProt:	$(findProtOBJS)
	$(CXX) -o $@ $(findProtOBJS) $(LIBS) 
	$(ECHO)

GeneMark2CDS.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
GeneMark2CDSOBJS=GeneMark2CDS.o $(SEQ_OBJ)
GeneMark2CDS:	$(GeneMark2CDSOBJS)
	$(CXX) -o $@ $(GeneMark2CDSOBJS) $(LIBS) 
	$(ECHO)

hmm_tc1.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DM_DIR)/numeric.hpp $(GEN_DIR)/hmm.hpp
hmm_tc1OBJS=hmm_tc1.o $(HMM_OBJ) 
hmm_tc1:	$(hmm_tc1OBJS)
	$(CXX) -o $@ $(hmm_tc1OBJS) $(LIBS) 
	$(ECHO)

hmm2prot.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp $(GEN_DIR)/hmm.hpp
hmm2protOBJS=hmm2prot.o $(HMM_OBJ) 
hmm2prot:	$(hmm2protOBJS)
	$(CXX) -o $@ $(hmm2protOBJS) $(LIBS) 
	$(ECHO)

hmmAddCutoff.o:  $(COMMON_HPP) 
hmmAddCutoffOBJS=hmmAddCutoff.o $(CPP_DIR)/common.o
hmmAddCutoff:	$(hmmAddCutoffOBJS)
	$(CXX) -o $@ $(hmmAddCutoffOBJS) $(LIBS) 
	$(ECHO)

hmmNAME2ACC.o:  $(COMMON_HPP) 
hmmNAME2ACCOBJS=hmmNAME2ACC.o $(CPP_DIR)/common.o
hmmNAME2ACC:	$(hmmNAME2ACCOBJS)
	$(CXX) -o $@ $(hmmNAME2ACCOBJS) $(LIBS) 
	$(ECHO)

hmmsearch2besthits.o:  $(COMMON_HPP) $(GEN_DIR)/hmm.hpp
hmmsearch2besthitsOBJS=hmmsearch2besthits.o $(HMM_OBJ) 
hmmsearch2besthits:	$(hmmsearch2besthitsOBJS)
	$(CXX) -o $@ $(hmmsearch2besthitsOBJS) $(LIBS) 
	$(ECHO)

kmerIndex_add.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
kmerIndex_addOBJS=kmerIndex_add.o $(SEQ_OBJ)
kmerIndex_add:	$(kmerIndex_addOBJS)
	$(CXX) -o $@ $(kmerIndex_addOBJS) $(LIBS) 
	$(ECHO)

kmerIndex_find.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
kmerIndex_findOBJS=kmerIndex_find.o $(SEQ_OBJ)
kmerIndex_find:	$(kmerIndex_findOBJS)
	$(CXX) -o $@ $(kmerIndex_findOBJS) $(LIBS) 
	$(ECHO)

kmerIndex_make.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
kmerIndex_makeOBJS=kmerIndex_make.o $(SEQ_OBJ)
kmerIndex_make:	$(kmerIndex_makeOBJS)
	$(CXX) -o $@ $(kmerIndex_makeOBJS) $(LIBS) 
	$(ECHO)

kmerIndex_stat.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
kmerIndex_statOBJS=kmerIndex_stat.o $(SEQ_OBJ)
kmerIndex_stat:	$(kmerIndex_statOBJS)
	$(CXX) -o $@ $(kmerIndex_statOBJS) $(LIBS) 
	$(ECHO)

mlst2hash.o:  $(COMMON_HPP) 
mlst2hashOBJS=mlst2hash.o $(CPP_DIR)/common.o
mlst2hash:	$(mlst2hashOBJS)
	$(CXX) -o $@ $(mlst2hashOBJS) $(LIBS) 
	$(ECHO)

mutation_tab.o:	$(COMMON_HPP) 
mutation_tabOBJS=mutation_tab.o $(CPP_DIR)/common.o
mutation_tab:	$(mutation_tabOBJS)
	$(CXX) -o $@ $(mutation_tabOBJS) $(LIBS) 
	$(ECHO)

orf2prot.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp
orf2protOBJS=orf2prot.o $(SEQ_OBJ)
orf2prot:	$(orf2protOBJS)
	$(CXX) -o $@ $(orf2protOBJS) $(LIBS) 
	$(ECHO)

prot_check.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot_checkOBJS=prot_check.o $(SEQ_OBJ)
prot_check:	$(prot_checkOBJS)
	$(CXX) -o $@ $(prot_checkOBJS) $(LIBS) 
	$(ECHO)

prot_fingerprint_grep.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot_fingerprint_grepOBJS=prot_fingerprint_grep.o $(SEQ_OBJ)
prot_fingerprint_grep:	$(prot_fingerprint_grepOBJS)
	$(CXX) -o $@ $(prot_fingerprint_grepOBJS) $(LIBS) 
	$(ECHO)

prot_grep_short.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot_grep_shortOBJS=prot_grep_short.o $(SEQ_OBJ)
prot_grep_short:	$(prot_grep_shortOBJS)
	$(CXX) -o $@ $(prot_grep_shortOBJS) $(LIBS) 
	$(ECHO)

prot2fingerprints.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot2fingerprintsOBJS=prot2fingerprints.o $(SEQ_OBJ)
prot2fingerprints:	$(prot2fingerprintsOBJS)
	$(CXX) -o $@ $(prot2fingerprintsOBJS) $(LIBS) 
	$(ECHO)

prot2triplets.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot2tripletsOBJS=prot2triplets.o $(SEQ_OBJ)
prot2triplets:	$(prot2tripletsOBJS)
	$(CXX) -o $@ $(prot2tripletsOBJS) $(LIBS) 
	$(ECHO)

prots2hmm_signature.o:  $(COMMON_HPP) $(GEN_DIR)/hmm.hpp
prots2hmm_signatureOBJS=prots2hmm_signature.o $(HMM_OBJ) 
prots2hmm_signature:	$(prots2hmm_signatureOBJS)
	$(CXX) -o $@ $(prots2hmm_signatureOBJS) $(LIBS) 
	$(ECHO)

prots_pair2stat.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DM_DIR)/numeric.hpp
prots_pair2statOBJS=prots_pair2stat.o $(SEQ_OBJ) $(DM_DIR)/numeric.o
prots_pair2stat:	$(prots_pair2statOBJS)
	$(CXX) -o $@ $(prots_pair2statOBJS) $(LIBS) 
	$(ECHO)

replaceFastaHeader.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
replaceFastaHeaderOBJS=replaceFastaHeader.o $(SEQ_OBJ)
replaceFastaHeader:	$(replaceFastaHeaderOBJS)
	$(CXX) -o $@ $(replaceFastaHeaderOBJS) $(LIBS) 
	$(ECHO)

replaceFastaNames.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
replaceFastaNamesOBJS=replaceFastaNames.o $(SEQ_OBJ)
replaceFastaNames:	$(replaceFastaNamesOBJS)
	$(CXX) -o $@ $(replaceFastaNamesOBJS) $(LIBS) 
	$(ECHO)

splitFasta.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
splitFastaOBJS=splitFasta.o $(SEQ_OBJ)
splitFasta:	$(splitFastaOBJS)
	$(CXX) -o $@ $(splitFastaOBJS) $(LIBS) 
	$(ECHO)

splitHmmLib.o:	$(COMMON_HPP) 
splitHmmLibOBJS=splitHmmLib.o $(CPP_DIR)/common.o
splitHmmLib:	$(splitHmmLibOBJS)
	$(CXX) -o $@ $(splitHmmLibOBJS) $(LIBS) 
	$(ECHO)

tblastn2orfs.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
tblastn2orfsOBJS=tblastn2orfs.o $(SEQ_OBJ)
tblastn2orfs:	$(tblastn2orfsOBJS)
	$(CXX) -o $@ $(tblastn2orfsOBJS) $(LIBS) 
	$(ECHO)

uniqSeq.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
uniqSeqOBJS=uniqSeq.o $(SEQ_OBJ)
uniqSeq:	$(uniqSeqOBJS)
	$(CXX) -o $@ $(uniqSeqOBJS) $(LIBS) 
	$(ECHO)

uniqProtRef.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
uniqProtRefOBJS=uniqProtRef.o $(SEQ_OBJ)
uniqProtRef:	$(uniqProtRefOBJS)
	$(CXX) -o $@ $(uniqProtRefOBJS) $(LIBS) 
	$(ECHO)

