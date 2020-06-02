#!/bin/csh -f

if ($# != 5) then
  echo "Blast and compute genome distances, save BLASTN output #1-#2.blast in #5/"
  echo "#1: genome 1 id"
  echo "#2: genome 2 id"
  echo "#3: Directory with DNA"
  echo "#4: megablast|tblastx"
  echo "#5: ANI|LANI (locus ANI)|cons"
  echo "Requires: #1 != #2; #3/{#1,#2}.{fa,n*}"
  exit 1
endif



# CUT, PROG
set CUT = "" 
switch ("$5")
  case "ANI":
    set CUT = "qseq sseq"  
    set PROG = "blast2ani 3000 -cut"  
    breaksw
  case "LANI":  # ??
    set PROG = "blast2ani 180"  
    breaksw
  case "cons":
    set L1 = `cat $3/$1.len`
    set L2 = `cat $3/$2.len`
    set PROG = "blast2cons 0 $L1 $L2"
    breaksw
  default:
    echo "Unknown distance $5"
    exit 1
    breaksw
endsw


if ("$4" == "megablast") then
  if (1) then
    blastn \
      -task megablast \
      -query $3/$1.fa \
      -db $3/$2.fa \
      -parse_deflines \
      -max_target_seqs 100000 \
      -xdrop_gap 150 \
      -xdrop_gap_final 150 \
      -penalty -1 \
      -gapopen 3 \
      -gapextend 1 \
      -dbsize 10000000 \
      -searchsp 10000000 \
      -outfmt "6 qseqid sseqid qstart qend sstart send nident length $CUT" \
      | $PROG  
      #> $1-$2.blastn
      # LANI: awk '$1 == $2'
    if ($?) exit 1
  else
    /panfs/pan1.be-md.ncbi.nlm.nih.gov/gpipe/ThirdParty/WUBLAST/2006-01-02/arch/i686/blastn $3/$2.fa $3/$1.fa \
      -gi \
      -B=100000 \
      -V=100000 \
      -gapw=1000 \
      -cpus=1 \
      -warnings -errors -notes \
      -Q 3 \
      -R 1 \
      -M 1 \
      -N -1 \
      -wordmask seg \
      -top \
      -mformat 2 \
      > $1-$2.wublastn 
    if ($?) exit 1
    cat $1-$2.wublastn | awk '{printf "%d %d %d %d %d %d %d %d\n", $1, $2, $18, $19, $21, $22, $8, $7;}' > $1-$2.blast
    if ($?) exit 1
  endif
else
  if (1) then
    /panfs/pan1.be-md.ncbi.nlm.nih.gov/gpipe/ThirdParty/WUBLAST/2006-01-02/arch/i686/tblastx $3/$1.fa $3/$2.fa \
      -gi \
      -B=100000 \
      -V=100000 \
      -gapw=1000 \
      -cpus=1 \
      -warnings -errors -notes \
      C=11 \
      -dbgcode 11 \
      -wordmask seg \
      -S 200 \
      -s2 200 \
      -W 5 \
      -mformat 2 \
      > tblastx/$1-$2
    if ($?) exit 1
  else
    tblastx \
      -query $3/$1.fa \
      -db $3/$2.fa \
      -parse_deflines \
      -max_target_seqs 100000 \
      -query_gencode 11 \
      -db_gencode 11 \
      -dbsize   1000000 \
      -searchsp 1000000 \
      -evalue 1e-5 \
      -word_size 5 \
      -outfmt '6 qseqid sseqid qstart qend sstart send nident length' \
      -out tblastx/$1-$2
    if ($?) exit 1
    #-query_gencode ??
    #-db_gencode ??
  endif
endif


