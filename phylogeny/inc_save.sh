#!/bin/bash --noprofile
source bash_common.sh
cd ~brovervv/code/cpp/phylogeny

set -x

distTree_inc_save.sh $PANFS/GenBank/bacteria/inc      inc/genome/bacteria
distTree_inc_save.sh $PANFS/GenBank/bacteria_type/inc inc/genome/bacteria_type
distTree_inc_save.sh $PANFS/GenBank/prok/inc          inc/genome/prokaryota
distTree_inc_save.sh $PANFS/GenBank/Chlorophyta/inc   inc/genome/Chlorophyta
distTree_inc_save.sh $PANFS/GenBank/Fungi/inc         inc/genome/Fungi
distTree_inc_save.sh $PANFS/GenBank/Protists/inc      inc/genome/Protists
distTree_inc_save.sh $PANFS/GenBank/Viridiplantae/inc inc/genome/Viridiplantae
distTree_inc_save.sh $PANFS/GenBank/Covid/inc         inc/virus/SARS-CoV-2
distTree_inc_save.sh $PANFS/RefSeq/Metazoa/inc        inc/genome/Metazoa-RefSeq
distTree_inc_save.sh $PANFS/RefSeq/Viridiplantae/inc  inc/genome/Viridiplantae-RefSeq
distTree_inc_save.sh $PANFS/marker/Fungi/18SrRNA/inc  inc/rRNA/Fungi/18SrRNA
distTree_inc_save.sh $PANFS/marker/Fungi/28SrRNA/inc  inc/rRNA/Fungi/28SrRNA
distTree_inc_save.sh $PANFS/marker/Fungi/ITS/inc      inc/rRNA/Fungi/ITS
distTree_inc_save.sh $PANFS/marker/SSU/inc            inc/rRNA/SSU
distTree_inc_save.sh $PANFS/marker/bacteria/inc       inc/rRNA/bacteria
distTree_inc_save.sh $PANFS/marker/5.8S/inc           inc/rRNA/5.8S


