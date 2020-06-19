create table List
(
  id numeric(20)  not null  primary key
);
go


create table ListC
(
  id varchar(256)  not null  primary key
);
go



create table Locus
(
  id varchar(12)  not null  constraint Locus_pk primary key

, dead bit  not null  default 0

, taxroot int  not null  
, gene varchar(16)  not null  

, in_tree bit  constraint Locus_in_tree_ck check (in_tree = 1)
, outlier varchar(16)  
, constraint Locus_outlier_in_tree_nl check (outlier is null or in_tree is null)
);
go
grant select on Locus to public;
create index Locus_gene_idx on Locus(gene);
create index Locus_outlier_idx on Locus(outlier);
go



create table Genome
(
  id int  not null  constraint Genome_pk primary key

, dead bit  not null  default 0

, total_len bigint  not null  constraint Genome_total_len_ck check (total_len >= 1000)
  
, tax_root int  not null

  -- Annotation
, annot_tried bit  not null  default 0
, constraint Genome_annot_tried_tax_root_ck check (annot_tried = 0 or tax_root is not null)

  -- Length >= 150 aa
, cdss int  constraint Genome_cdss_ck check (cdss >= 0)
    -- cdss is null and annot_tried = 1 <=> annotation failed
    -- from GenomeHash => only long sequences
, constraint Genome_annot_tried_cdss_nl check (annot_tried = 1 or cdss is null)
, prots int  constraint Genome_prots_ck check (prots >= 0)
    -- from GenomeHash => only long sequences
    -- annot_tried = 1 and prots is null: depends on N50
, constraint Genome_annot_tried_prots_ck check (annot_tried = 1 or prots is null)
, constraint Genome_prots_nl check (prots is null or cdss is not null)
, constraint Genome_cdss_prots_ck check (cdss >= prots)
, gene_length_ave  as total_len * 1.0 / cdss
    -- Average long gene length

, outlier varchar(64)   

  -- DNA statistics
, long_len int    constraint Genome_long_len check (long_len > 0) 
    -- DNA sequences of length >= 1000 bp
, constraint Genome_long_len_nl check (long_len is null or prots is not null)
, constraint Genome_total_long_len_ck check (long_len <= total_len)
, ambig real   constraint Genome_ambig_ck check (ambig between 0 and 1)
, constraint Genome_ambig_nl check (long_len is null and ambig is null or long_len is not null and ambig is not null)
, gc real   constraint Genome_gc_ck check (gc between 0 and 1)
, constraint Genome_gc_nl check (long_len is null and gc is null or long_len is not null and gc is not null)
, repeat8 real  constraint Genome_repeat8_ck  check (repeat8 between 0 and 1)
, constraint Genome_repeat8_nl check (long_len is null and repeat8 is null or long_len is not null and repeat8 is not null)

  -- Tree
, in_tree bit  constraint Genome_in_tree_ck check (in_tree = 1)
    -- 1:    in tree
    -- null: in new/
, constraint Genome_prots_in_tree_ck check (prots is not null or in_tree is null)
, genogroup int  constraint Genome_genogroup_fk references Genome(id)
    -- Genospecies or genosubspecies
, constraint Genome_genogroup_prots_nl check (genogroup is null or prots is not null)
);
go
grant select on Genome to public;
create index Genome_genogroup_idx on Genome(genogroup);
go


create function normal_distr (@x real, @mean real, @sd real) returns real
as
begin
  declare @y real = (@x - @mean) / (@sd * sqrt (2));
  declare @sign real = 1;
  if @y < 0
  begin
    set @y = - @y;
    set @sign = -1;
  end;
  -- https://en.wikipedia.org/wiki/Error_function
  -- Approximate
  declare @sqrt_pi real = sqrt (pi());
  declare @ey real = exp (- power (@y,2));
  declare @erf real = 2 / @sqrt_pi * sqrt (1 - @ey) * (@sqrt_pi / 2 + 31.0/200 * @ey - 341.0/8000 * power(@ey, 2));
  -- https://en.wikipedia.org/wiki/Normal_distribution
  return 0.5 * (1 + @sign * @erf);
end;
go


create function quantil2std_normal (@quantil real) returns real
as
begin
  if @quantil < 0
    return 0.0 / 0.0;  -- error
  declare @x_lo real = -1000;  -- PAR
  declare @x_hi real =  1000;  -- PAR
  declare @x real;
  declare @d real;
  -- Binary search
  while @x_hi - @x_lo > 0.001  -- PAR
  begin
    set @x = (@x_lo + @x_hi) / 2;
    set @d = dbo.normal_distr (@x, 0, 1);
    if @d > @quantil
      set @x_hi = @x;
    else
      set @x_lo = @x;
  end;
  return @x;
end;
go


create procedure Genogroup2outliers (@tax_root int, @use_weak bit)
-- Input: ListC: <Genome.id>\t<genogroup>
-- Print: List of outlier Genome.id's 
-- Converges
as
begin
  set nocount on;

  update Genome
    set genogroup = null
    where tax_root = @tax_root;
  update Genome
    set genogroup = substring(ListC.id,charindex(char(9),ListC.id)+1,len(ListC.id))
    from      ListC
         join Genome on Genome.id = substring(ListC.id,1,charindex(char(9),ListC.id)-1)
    where Genome.tax_root = @tax_root;  -- redundant

--drop table #Genogroup;
  select genogroup id
       , count(*) c
       , avg(total_len * 1.0) total_len_avg
       , stdev(total_len * 1.0) total_len_sd
       , avg(gc) gc_avg
       , stdev(gc) gc_sd
       , avg(repeat8) repeat8_avg
       , stdev(repeat8) repeat8_sd
       , avg(gene_length_ave) gene_length_ave_avg
       , stdev(gene_length_ave) gene_length_ave_sd
    into #Genogroup
    from Genome
    where     tax_root = @tax_root
          and genogroup is not null
          and outlier is null
          and in_tree = 1
    group by genogroup
    having count(*) >= 10;  -- PAR
  alter table #Genogroup alter column id int not null;
  alter table #Genogroup add primary key (id);

  declare @evalue real = 0.01;  -- PAR
  alter table #Genogroup add threshold real  check (threshold > 0);
  update #Genogroup
    set threshold = - dbo.quantil2std_normal (0.5 * @evalue / c);  -- two-tail test

  -- Genome.outlier
  update Genome
    set outlier = 'gene_length_ave: ' + case when Genome.gene_length_ave > #Genogroup.gene_length_ave_avg then 'large' else 'small' end 
    from      Genome
         join #Genogroup on #Genogroup.id = Genome.genogroup
    where     Genome.tax_root = @tax_root
          and Genome.outlier is null
          and in_tree = 1
          and abs(Genome.gene_length_ave - #Genogroup.gene_length_ave_avg) > #Genogroup.threshold * #Genogroup.gene_length_ave_sd;
  update Genome
    set outlier = 'repeat8: ' + case when Genome.repeat8 > #Genogroup.repeat8_avg then 'large' else 'small' end 
    from      Genome
         join #Genogroup on #Genogroup.id = Genome.genogroup
    where     Genome.tax_root = @tax_root
          and Genome.outlier is null
          and in_tree = 1
          and abs(Genome.repeat8 - #Genogroup.repeat8_avg) > #Genogroup.threshold * #Genogroup.repeat8_sd;
  if @use_weak = 1
  begin
    update Genome
      set outlier = 'total_len: ' + case when Genome.total_len > #Genogroup.total_len_avg then 'large' else 'small' end 
      from      Genome
           join #Genogroup on #Genogroup.id = Genome.genogroup
      where     Genome.tax_root = @tax_root
            and Genome.outlier is null
            and in_tree = 1
            and abs(Genome.total_len - #Genogroup.total_len_avg) > #Genogroup.threshold * #Genogroup.total_len_sd;
    update Genome
      set outlier = 'gc: ' + case when Genome.gc > #Genogroup.gc_avg then 'large' else 'small' end 
      from      Genome
           join #Genogroup on #Genogroup.id = Genome.genogroup
      where     Genome.tax_root = @tax_root
            and Genome.outlier is null
            and in_tree = 1
            and abs(Genome.gc - #Genogroup.gc_avg) > #Genogroup.threshold * #Genogroup.gc_sd;
  end;

  select id
    from Genome
    where     tax_root = @tax_root
          and outlier is not null
          and in_tree = 1;

  set nocount off;
end;
go


create table GenomeHash
(
  [hash] numeric(20)  not null
, [type] char(3)  not null  
, genome int  not null  
);
create unique index GenomeHash_hash_type_uq on GenomeHash([hash], [type], genome);
go


create table FreqHash
(
  tax_root int  not null  
, [type] char(3)  not null  
, [hash] numeric(20)  not null
, primary key ([hash], tax_root, [type])
);
go
