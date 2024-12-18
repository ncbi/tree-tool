#!/bin/bash --noprofile
THIS=$( realpath $( dirname $0 ) )
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Load NCBI taxonomy data into the database #2 on the Microsoft SQL Server server #1, schema #3"
  echo "Output tables: Taxrank, Tax"
  echo 'This command must connect to the database without requesting a password and permit DDL operations "sqsh-ms -S $SERVER -D $DATABASE"'
  echo "#1: Server name"
  echo "#2: Database name"
  echo "#3: Schema name"
  echo "Time: 5 min."
  exit 1
fi
SERVER=$1
DB=$2
SCHEMA=$3


#if false; then
TMP=$( mktemp )
comment $TMP


mkdir $TMP.dmp
cd $TMP.dmp


section "Downloading data"
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz  -O $TMP.dmp/taxdump.tar.gz
#    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz ??
tar -xzf taxdump.tar.gz

if [ ! -e nodes.dmp ]; then
  error "No nodes.dmp"
fi

if [ ! -e names.dmp ]; then
  error "No names.dmp"
fi


sqsh-ms  -S $SERVER  -D $DB << EOT > $TMP.schema
  select 1 
    from sys.schemas 
    where [name] = '$SCHEMA'
  go -m bcp
EOT
if [ ! -s $TMP.schema ]; then
  sqsh-ms  -S $SERVER  -D $DB << EOT 
    create schema $SCHEMA
    go
EOT
fi


section "Creating tables"
sqsh-ms  -S $SERVER  -D $DB << EOT 
  begin try
    create table $SCHEMA.Taxrank
    (
      id varchar(16)  not null  primary key
    , num int  not null  check (num >= -1)
    , sub_rank_max tinyint  not null  default 0 
    )
    insert into $SCHEMA.Taxrank (num, id) values
      (-1, 'all'),
      ( 0, 'superkingdom'),
      ( 1, 'kingdom'),
      ( 2, 'subkingdom'),
      ( 3, 'superphylum'),
      ( 4, 'phylum'),
      ( 5, 'subphylum'),
      ( 6, 'superclass'),
      ( 7, 'class'),
      ( 8, 'subclass'),
      ( 9, 'infraclass'),
      (10, 'cohort'),
      (11, 'subcohort'),
      (12, 'superorder'),
      (13, 'order'),
      (14, 'suborder'),
      (15, 'infraorder'),
      (16, 'parvorder'),
      (17, 'superfamily'),
      (18, 'family'),
      (19, 'subfamily'),
      (20, 'tribe'),
      (21, 'subtribe'),
      (22, 'genus'),
      (23, 'subgenus'),
      (24, 'species group'),
      (25, 'species subgroup'),
      (26, 'species'),
      (27, 'subspecies'),
      (28, 'varietas'),
      (29, 'forma')
      -- 'section', 'series'
    grant select on $SCHEMA.Taxrank to public
    
    create table $SCHEMA.Tax
    (
      id int  not null  primary key  check (id > 0)  
    , parent int  
    , rank_name varchar(16)
    , [rank] int
    , [name] varchar(256)  not null  
         check (    not [name] like ' %'
                and not [name] like '% '
                and not [name] like '%;%'
               )
        -- unique scientific name
    , revised as  case when substring([name],1,1) = '[' then 1 else 0 end
    , sub_rank tinyint  not null  default 0  
    , endosymbiont bit  not null  default 0
    , strain bit  not null  default 0
    , unclassified bit  not null  default 0
    , depth int  check (depth >= 0)
    )
    grant select on $SCHEMA.Tax to public
  end try
  begin catch
    print error_message()
    raiserror ('error', 11, 1)
  end catch
  go
EOT


section "Populating Tax"
awk -F '\t' '$7 =="scientific name"' names.dmp | awk -F '\t' '$5 != ""' | cut -f 1,5 >  $TMP.names
awk -F '\t' '$7 =="scientific name"' names.dmp | awk -F '\t' '$5 == ""' | cut -f 1,3 >> $TMP.names
set +o errexit
grep '~' $TMP.names  # needed for Tax1000.sh
S=$?
set -o errexit
if [ $S != 1 ]; then
  exit 1
fi
sed "s/'/''/g" $TMP.names | sort > $TMP.names-sorted

cut -f 1,3,5 nodes.dmp | sort > $TMP.nodes
join  -1 1  -2 1  -t $'\t'  $TMP.nodes $TMP.names-sorted > $TMP.taxdata

mkdir $TMP.Tax
$THIS/../splitList $TMP.taxdata 1000 $TMP.Tax
mkdir $TMP.sql
$THIS/../trav $TMP.Tax  -threads 10  -step 1  "$THIS/Tax1000.sh $SCHEMA %d/%f > $TMP.sql/%f"
mkdir $TMP.out
$THIS/../trav $TMP.sql  -threads 10  -step 1  "sqsh-ms -S $SERVER  -D $DB  -i %d/%f > $TMP.out/%f"
set +o errexit
grep -v '^([1-9].* rows affected)$' $TMP.out/*
S=$?
set -o errexit
if [ $S != 1 ]; then
  exit 1
fi


section "Indexing"
sqsh-ms  -S $SERVER  -D $DB << EOT 
  begin try
    create index Tax_parent_idx on $SCHEMA.Tax(parent)
    update $SCHEMA.Tax 
      set parent = null
      where id = 1
    alter table $SCHEMA.Tax add constraint Tax_fk foreign key (parent) references $SCHEMA.Tax(id)
    update $SCHEMA.Tax
      set rank_name = null
      where not rank_name in (select Taxrank.id from $SCHEMA.Taxrank)
    alter table $SCHEMA.Tax add constraint Tax_rank_fk foreign key (rank_name) references $SCHEMA.Taxrank(id)
    create unique index Tax_uq on $SCHEMA.Tax([name])
    create index Tax_depth_idx on $SCHEMA.Tax(depth)
    update Tax
      set [rank] = Taxrank.num
      from      $SCHEMA.Tax
           join $SCHEMA.Taxrank on Taxrank.id = Tax.rank_name
    create index Tax_rank_idx on $SCHEMA.Tax([rank])
  end try
  begin catch
    print error_message()
    raiserror ('error', 11, 1)
  end catch
  go
EOT


section "Creating stored functions"
sqsh-ms  -S $SERVER  -D $DB << EOT 
  create function $SCHEMA.tax_rank_sub_rank2tax (@tax int, @rank int, @sub_rank tinyint) returns int
  -- Traverse Tax lineage from @tax up to (@rank,@sub_rank)
  -- Return: Tax.id, null (root) or -1 (unclassified)
  as 
  begin   
    declare @parent int  
    declare @found bit
    declare @unclassified bit = 0
    while 1=1
    begin
      set @found = null
      select @parent = parent
           , @found = case when [rank] = @rank and sub_rank = @sub_rank and strain = 0 then 1 else 0 end
           , @unclassified = unclassified
        from $SCHEMA.Tax
        where id = @tax
      if @found is null
        return null 
      if @found = 1 
      begin
        if @unclassified = 1
          return -1
        return @tax
      end
      if @parent is null
        return null
      set @tax = @parent
    end
    return null
  end
  go


  create function $SCHEMA.tax_rank2tax (@tax int, @rank int) returns int
  as
  begin
    return $SCHEMA.tax_rank_sub_rank2tax (@tax, @rank, 0)
  end
  go


  create function $SCHEMA.tax_rank_sub_rank2name (@tax int, @rank int, @sub_rank tinyint) returns varchar(256)
  -- Return: null (root), '_OTHER_' (unclassified) or name
  as
  begin
    declare @ancestor_tax int
    select @ancestor_tax = $SCHEMA.tax_rank_sub_rank2tax (@tax, @rank, @sub_rank)
    if @ancestor_tax is null
      return null
    if @ancestor_tax = -1
      return '_OTHER_'
    declare @name varchar(256)
    select @name = [name]
      from $SCHEMA.Tax
      where id = @ancestor_tax
    return @name
  end
  go


  create procedure $SCHEMA.tax2phen (@tax int, @rank_min int)
  -- Print: <NN>-<NN>:{<tax_name>|_OTHER_}
  as
  begin
    declare @rank int
    select @rank = Tax.[rank]
      from $SCHEMA.Tax
      where Tax.id = @tax
    if @rank is null
      raiserror ('No record in $SCHEMA.Tax', 11, @tax)

    declare RankCur cursor local for
      select num, sub_rank_max
        from $SCHEMA.Taxrank
        where num between @rank_min and @rank  -- PAR
        order by num
    open RankCur
    declare @taxrank int
    declare @sub_rank_max tinyint
    declare @taxname varchar(256)
    declare @taxname_old varchar(128) = ''
    declare @sub_rank tinyint
    declare @name_prefix varchar(128)
    while 1=1
    begin
      fetch next from RankCur into @taxrank, @sub_rank_max
      if @@fetch_status != 0
        break
      set @sub_rank = 0
      while @sub_rank <= @sub_rank_max
      begin
        if @taxrank >= 100
          raiserror ('@taxrank >= 100', 11, @tax)
        if @sub_rank >= 100
          raiserror ('@sub_rank >= 100', 11, @tax)
        set @taxname = $SCHEMA.tax_rank_sub_rank2name (@tax, @taxrank, @sub_rank)
        set @name_prefix = replace(str(@taxrank, 2),' ','0') + '-' + replace(str(@sub_rank,2),' ','0') + ':'
        if @taxname is not null
        begin
          print @name_prefix + @taxname
          set @taxname_old = @taxname
        end
        else
        begin
          if @taxname_old != '' 
            print @name_prefix + @taxname_old  
        end
        set @sub_rank = @sub_rank + 1
      end
    end
    close RankCur
    deallocate RankCur
  end
  go
  grant execute on $SCHEMA.tax2phen to public
  go
EOT


section "Finishing"
sqsh-ms  -S $SERVER  -D $DB << EOT 
  begin try
    set nocount on

    -- $SCHEMA.Tax.{rank_name,[rank],sub_rank}
    update $SCHEMA.Tax
      set rank_name = 'all'
        , [rank] = -1
      where     parent is null
            and rank_name is null
    while 1=1
    begin
      update Tax
        set rank_name = P.rank_name
          , [rank]    = P.[rank]
          , sub_rank  = P.sub_rank + 1
        from      $SCHEMA.Tax
             join $SCHEMA.Tax P on P.id = Tax.parent
        where     Tax.rank_name is null
              and P.rank_name is not null
      if @@rowcount = 0
        break
    end
    update $SCHEMA.Tax
      set rank_name = 'all'
        , [rank] = -1
      where rank_name is null

    -- Tax.depth
    update $SCHEMA.Tax
      set depth = 0
      where parent is null
    while 1=1
    begin
      update Tax
        set depth = P.depth + 1
        from      $SCHEMA.Tax
             join $SCHEMA.Tax P on P.id = Tax.parent
        where     Tax.depth is null
              and P.depth is not null
      if @@rowcount = 0
        break
    end

    -- Taxrank.sub_rank_max
  	update Taxrank
  	  set sub_rank_max = T.sub_rank_max
  	  from      $SCHEMA.Taxrank
  	       join (select rank_name, max(sub_rank) sub_rank_max
                   from $SCHEMA.Tax
                   group by rank_name
                ) T on T.rank_name = Taxrank.id

    -- Tax.endosymbiont
    update $SCHEMA.Tax
      set endosymbiont = 1
      where id in
              ( 9       -- Buchnera aphidicola
              , 186490  -- Candidatus Baumannia cicadellinicola
              , 336810  -- Candidatus Sulcia muelleri
              , 813     -- Chlamydia trachomatis
              , 780     -- Rickettsia
              , 138074  -- Serratia symbiotica
              , 51229	  -- Wigglesworthia glossinidia
              )
    update $SCHEMA.Tax
      set endosymbiont = 1
      where lower(name) like '%endosymbiont%'
    update $SCHEMA.Tax
      set endosymbiont = 1
      where lower(name) like '%phytoplasma%'
    while 1=1
    begin
      update Tax
        set endosymbiont = 1
        from      $SCHEMA.Tax
             join $SCHEMA.Tax P on P.id = Tax.parent
        where     Tax.endosymbiont = 0
              and P.endosymbiont = 1
      if @@rowcount = 0
        break
    end

    -- Tax.unclassified
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '%unclassified%'
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '%unassigned%'
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '%unidentified%'
    update $SCHEMA.Tax
      set unclassified = 1
      where [name] like '% candidate phyla'
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '%incertae sedis%'
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '%symbionts'
    update $SCHEMA.Tax
      set unclassified = 1
      where     name like '% sp.%'
            and not [name] like '% f. sp.%'  -- "Forma speciales"
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '% taxa'
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '% taxon %'
    update $SCHEMA.Tax
      set unclassified = 1
      where lower(name) like '%uncultured%'
    update $SCHEMA.Tax
      set unclassified = 1
      where name like 'bacterium%'
    update $SCHEMA.Tax
      set unclassified = 1
      where name like '% bacterium%'
    update $SCHEMA.Tax
      set unclassified = 1
      where name like 'actinobacterium%'
    update $SCHEMA.Tax
      set unclassified = 1
      where name like '% actinobacterium%'
    update $SCHEMA.Tax
      set unclassified = 1
      where name like 'proteobacterium%'
    update $SCHEMA.Tax
      set unclassified = 1
      where name like '% proteobacterium%'
    update $SCHEMA.Tax
      set unclassified = 1
      where     revised = 0
            and sub_rank = 0
            and endosymbiont = 0
            and not [name] like 'candidate %'
            and not [name] like '% group'
            and not [name] like '''%'
            and substring([name],1,1) = lower(substring([name],1,1))
    update $SCHEMA.Tax
      set unclassified = 1
      where     id = 46583  -- Candida tanzawaensis
            and $SCHEMA.tax_rank2tax (id, 22/*genus*/) is null


    -- Tax.strain
    update $SCHEMA.Tax
      set strain = 1
      where     [rank] >= 26
            and replace([name],'Candidatus ','') like '% % %'
            and (   [name] like '% str. %'
                 or [name] like '%[ (]strain %'
                 or [name] like '%[ (]NRCC%'
                 or [name] like '%[ (]ATCC%'
                 or [name] like '%[ (]PCC%'
                 or [name] like '%[ (]DSM%'
                 or [name] like '%[ (]NRRL%'
                 or [name] like '%[ (]CDC%'
                 or [name] like '%[ (]UBA%'
                 or [name] like '%[ (]SCGC%'
                 or [name] like '%[ (]JSM%'
                 or [name] like '%[ (]JCM%'
                 or [name] like '%[ (]IFO%'
                 or [name] like '%[ (]NBRC%'
                )
    update $SCHEMA.Tax
      set strain = 1
      where     [rank] >= 26
            and replace([name],'Candidatus ','') like '% % %'
            and id in (103721, 1183427, 1300253, 1004786, 1300255, 1104605, 2086579, 864596, 1384461,
                       1385724, 1379303, 2043164, 2043165, 2043166, 2043160, 1702170, 212717, 688245, 1204414, 2043169, 1380412, 1379303,
                       679895, 595496, 668369, 718252, 298653, 272568, 657308, 1255570, 478009, 1354721, 1055530, 1755239, 1706231, 1962180, 
                       272620, 378753, 1263082, 393127, 637915, 1201294, 876269, 1807683, 1807684, 1807685, 1807686, 1807687, 1807688, 1807689,
                       1960879, 1727164, 1158338, 1158345, 754252, 281093, 1078464, 1749078, 1798213, 1165861, 1125979, 234621, 193079, 1298918,
                       2043162, 1112349, 1297534, 869303, 869304, 869306, 869307, 1882757, 1984801, 2043167, 2043168, 1313298, 1629664, 223789,
                       679716, 223926, 223926)
    while 1=1
    begin
      update C
        set strain = 1
        from      $SCHEMA.Tax P
             join $SCHEMA.Tax C on C.parent = P.id
        where     P.strain = 1
              and C.strain = 0
      if @@rowcount = 0
        break
    end


    -- Integrity
    declare @n int
    
    select @n = count(*)
      from $SCHEMA.Tax
      where     parent is null
            and id != 1
    if @n > 0
      raiserror ('$SCHEMA.Tax.parent', 11, @n)

    select @n = count(*)
      from $SCHEMA.Tax
      where    [rank] is null
            or rank_name is null
    if @n > 0
      raiserror ('$SCHEMA.Tax.rank/rank_name is null', 11, @n)

    select @n = count(*)
      from $SCHEMA.Tax
      where depth is null
    if @n > 0
      raiserror ('Tax.depth', 11, @n)

    select @n = count(*)
      from $SCHEMA.Tax
      where not (strain = 0 or [rank] >= 26)
    if @n > 0
      raiserror ('Tax.strain rank is too high', 11, @n)

    select @n = count(*)
      from      $SCHEMA.Tax C
           join $SCHEMA.Tax P on P.id = C.parent
      where     C.[rank] > P.[rank]
            and C.sub_rank >= 1
    if @n > 0
      raiserror ('Tax/Parent.[rank]', 11, @n)
  end try
  begin catch
    print error_message()
    raiserror ('error', 11, 1)
  end catch
  go -m bcp 
EOT


rm -r $TMP*


success
