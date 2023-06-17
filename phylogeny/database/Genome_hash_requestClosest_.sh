#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../../bash_common.sh
if [ $# != 7 ]; then
  echo "Print Genome.id's approximately closest to #1, where Genome.in_tree = 1"
  echo "#1: SQL Server name"
  echo "#2: Database with tables Genome, GenomeHash, FreqHash"
  echo "#3: path in Universal Naming Convention to the bulk directory"
  echo "#4: Genome.id in #2"
  echo "#5: Taxroot.id"
  echo "#6: Hash type"
  echo "#7: file with hashes in the bulk directrory"
  exit 1
fi
SERVER=$1
DATABASE=$2
BULK_REMOTE=$3
GENOME=$4
TAXROOT=$5
TYPE=$6
FILE=$7


#set -x


sqsh-ms  -S $SERVER  -D $DATABASE << EOT | sed 's/|$//1' | grep -vx $GENOME
  -- Time: O(h log (h n) + sort(h k)), where h - # hashes per genome, k - # genomes per hash
  --       k = O(1) artificially
  -- Time to populate GenomeHash: O(sort(h n))+
  begin
    create table #H ([hash] numeric(20)  not null)
    bulk insert #H from '$BULK_REMOTE\\$FILE'
      with
      (
        BATCHSIZE = 100000
      , CHECK_CONSTRAINTS
      , FIRE_TRIGGERS 
      , fieldterminator = '|'
      )

    select #H.[hash]  
      into #A
      from           #H
           left join FreqHash (nolock) on     FreqHash.[hash] = #H.[hash]
                                          and FreqHash.[type] = '$TYPE'
                                          and FreqHash.tax_root = $TAXROOT
      where FreqHash.[hash] is null

    select GenomeHash.genome, count(*) c
      into #G
      from      #A
           join GenomeHash (nolock) on     GenomeHash.[hash] = #A.[hash]
                                       and GenomeHash.[type] = '$TYPE'
      group by GenomeHash.genome

    select top 100/*PAR*/ #G.genome
      from      #G
           join Genome (nolock) on Genome.id = #G.genome
      where     Genome.in_tree = 1
          --and #G.genome != $GENOME
            and Genome.tax_root = $TAXROOT
      order by #G.c desc, #G.genome /*tie resolution*/
  end
  go -m bcp
EOT


