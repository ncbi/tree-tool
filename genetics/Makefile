CODE_DIR=$(shell pwd)/../..
include ../MakeRules
CPPFLAGS=$(OptCPPFLAGS)
#CPPFLAGS=$(DebugCPPFLAGS)

ifdef AT_NCBI
  CPPFLAGS += -Werror
endif



############################### Programs ################################

all:	\
  align2html \
  asm_gap \
  asnt \
	blastmat \
	blast_best_hits \
	blastp_merge \
	contig2read_coverage \
	dna2orfs \
	dna2stat \
	dna_complexity \
  dna_consensus \
	dna_coverage \
	dna_cut \
	dna_find \
	dna_gc_skew \
	dna_rev \
	dna_trim \
	extractFastaDna \
	extractFastaProt \
	fasta_prefix \
	fasta2GC \
	fasta2hash \
	fasta2len \
	fastaAddGi \
	fastaDna2pairs \
	genbank_grep \
	GeneMark2CDS \
	hmm2prot \
	hmm_tc1 \
	hmmAddCutoff \
	hmmExtract \
	hmmNAME2ACC \
	hmmsearch2besthits \
	hmmSplit \
	interSeq \
	islander \
	kmerIndex_add \
	kmerIndex_find \
	kmerIndex_make \
	kmerIndex_stat \
	mlst2hash \
	mutation_dna2prot \
	mutation_tab \
	orf2prot \
	prot2fingerprints \
	prot2triplets \
	prot_check \
	prot_consensus \
	prot_find \
	prot_fingerprint_grep \
	prot_grep_short \
	prots2hmm_signature \
	prots_pair2stat \
	replaceFastaHeader \
	replaceFastaNames \
	seq_print \
	splitFasta \
	tblastn2frameshift \
	tblastn2orfs \
	uniqProtRef \
	uniqSeq


align2html.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DM_HPP) 
align2htmlOBJS=align2html.o $(GEN_DIR)/seq.o $(DM_OBJ)
align2html:	$(align2htmlOBJS)
	$(CXX) -o $@ $(align2htmlOBJS) $(LIBS) 
	$(ECHO)

asm_gap.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp $(DM_HPP) 
asm_gapOBJS=asm_gap.o $(GEN_DIR)/seq.o $(DM_OBJ)
asm_gap:	$(asm_gapOBJS)
	$(CXX) -o $@ $(asm_gapOBJS) $(LIBS) 
	$(ECHO)

asnt.o:  $(COMMON_HPP) 
asntOBJS=asnt.o $(CPP_DIR)/common.o
asnt:	$(asntOBJS)
	$(CXX) -o $@ $(asntOBJS) $(LIBS) 
	$(ECHO)

blastmat.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
blastmatOBJS=blastmat.o $(SEQ_OBJ)
blastmat:	$(blastmatOBJS)
	$(CXX) -o $@ $(blastmatOBJS) $(LIBS) 
	$(ECHO)

blast_best_hits.o:  $(COMMON_HPP) 
blast_best_hitsOBJS=blast_best_hits.o $(CPP_DIR)/common.o
blast_best_hits:	$(blast_best_hitsOBJS)
	$(CXX) -o $@ $(blast_best_hitsOBJS) $(LIBS) 
	$(ECHO)

blastp_merge.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp $(DM_DIR)/matrix.hpp
blastp_mergeOBJS=blastp_merge.o $(CPP_DIR)/common.o $(DM_DIR)/numeric.o $(DM_DIR)/matrix.o
blastp_merge:	$(blastp_mergeOBJS)
	$(CXX) -o $@ $(blastp_mergeOBJS) $(LIBS) 
	$(ECHO)

contig2read_coverage.o:  $(COMMON_HPP) $(DM_DIR)/numeric.hpp $(TSV_DIR)/tsv.hpp
contig2read_coverageOBJS=contig2read_coverage.o $(CPP_DIR)/common.o $(DM_DIR)/numeric.o $(TSV_DIR)/tsv.o
contig2read_coverage:	$(contig2read_coverageOBJS)
	$(CXX) -o $@ $(contig2read_coverageOBJS) $(LIBS) 
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

dna_complexity.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_complexityOBJS=dna_complexity.o $(SEQ_OBJ)
dna_complexity:	$(dna_complexityOBJS)
	$(CXX) -o $@ $(dna_complexityOBJS) $(LIBS) 
	$(ECHO)

dna_consensus.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_consensusOBJS=dna_consensus.o $(SEQ_OBJ)
dna_consensus:	$(dna_consensusOBJS)
	$(CXX) -o $@ $(dna_consensusOBJS) $(LIBS) 
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

dna_find.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_findOBJS=dna_find.o $(SEQ_OBJ)
dna_find:	$(dna_findOBJS)
	$(CXX) -o $@ $(dna_findOBJS) $(LIBS) 
	$(ECHO)

dna_gc_skew.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
dna_gc_skewOBJS=dna_gc_skew.o $(SEQ_OBJ)
dna_gc_skew:	$(dna_gc_skewOBJS)
	$(CXX) -o $@ $(dna_gc_skewOBJS) $(LIBS) 
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

fasta_prefix.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
fasta_prefixOBJS=fasta_prefix.o $(SEQ_OBJ)
fasta_prefix:	$(fasta_prefixOBJS)
	$(CXX) -o $@ $(fasta_prefixOBJS) $(LIBS) 
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

genbank_grep.o:	$(COMMON_HPP) $(GEN_DIR)/genbank_text.hpp
genbank_grepOBJS=genbank_grep.o $(CPP_DIR)/common.o $(GEN_DIR)/genbank_text.o
genbank_grep:	$(genbank_grepOBJS)
	$(CXX) -o $@ $(genbank_grepOBJS) $(LIBS) 
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

hmmExtract.o:	$(COMMON_HPP) 
hmmExtractOBJS=hmmExtract.o $(CPP_DIR)/common.o
hmmExtract:	$(hmmExtractOBJS)
	$(CXX) -o $@ $(hmmExtractOBJS) $(LIBS) 
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

hmmSplit.o:	$(COMMON_HPP) 
hmmSplitOBJS=hmmSplit.o $(CPP_DIR)/common.o
hmmSplit:	$(hmmSplitOBJS)
	$(CXX) -o $@ $(hmmSplitOBJS) $(LIBS) 
	$(ECHO)

interSeq.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
interSeqOBJS=interSeq.o $(SEQ_OBJ)
interSeq:	$(interSeqOBJS)
	$(CXX) -o $@ $(interSeqOBJS) $(LIBS) 
	$(ECHO)

islander.o:  $(COMMON_HPP) $(TSV_DIR)/tsv.hpp
islanderOBJS=islander.o $(CPP_DIR)/common.o $(TSV_DIR)/tsv.o
islander:	$(islanderOBJS)
	$(CXX) -o $@ $(islanderOBJS) $(LIBS) 
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

mutation_dna2prot.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp
mutation_dna2protOBJS=mutation_dna2prot.o $(SEQ_OBJ)
mutation_dna2prot:	$(mutation_dna2protOBJS)
	$(CXX) -o $@ $(mutation_dna2protOBJS) $(LIBS) 
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

prot_consensus.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot_consensusOBJS=prot_consensus.o $(SEQ_OBJ)
prot_consensus:	$(prot_consensusOBJS)
	$(CXX) -o $@ $(prot_consensusOBJS) $(LIBS) 
	$(ECHO)

prot_find.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
prot_findOBJS=prot_find.o $(SEQ_OBJ)
prot_find:	$(prot_findOBJS)
	$(CXX) -o $@ $(prot_findOBJS) $(LIBS) 
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

seq_print.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
seq_printOBJS=seq_print.o $(SEQ_OBJ)
seq_print:	$(seq_printOBJS)
	$(CXX) -o $@ $(seq_printOBJS) $(LIBS) 
	$(ECHO)

splitFasta.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
splitFastaOBJS=splitFasta.o $(SEQ_OBJ)
splitFasta:	$(splitFastaOBJS)
	$(CXX) -o $@ $(splitFastaOBJS) $(LIBS) 
	$(ECHO)

tblastn2frameshift.o:  $(COMMON_HPP) 
tblastn2frameshiftOBJS=tblastn2frameshift.o $(CPP_DIR)/common.o
tblastn2frameshift:	$(tblastn2frameshiftOBJS)
	$(CXX) -o $@ $(tblastn2frameshiftOBJS) $(LIBS) 
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

virus2host.o:  $(COMMON_HPP) $(GEN_DIR)/seq.hpp 
virus2hostOBJS=virus2host.o $(SEQ_OBJ)
virus2host:	$(virus2hostOBJS)
	$(CXX) -o $@ $(virus2hostOBJS) $(LIBS) 
	$(ECHO)

