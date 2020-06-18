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



