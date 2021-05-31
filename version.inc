#define VERSION "1.9.3"


// 1.9.3  05/30/2021  align.cpp: ambiguous aa J,U,O are replaced by X, selfScore is using X
//                    -prefer_gain is removed
// 1.9.2  05/18/2021  optimizeReinsert() invokes optimizeLargeSubgraphs(&changes)
// 1.9.1  05/15/2021  inc/threads
//                    distTree_inc_new.sh -reinsert does not process inc/new/
// 1.8.1  05/11/2021  inc/dissim_full
// 1.7.3  04/23/2021  combine_dissims.cpp: weighted combination of all dissimilarities
//                    str2real() bug in reading nan's
// 1.7.2  04/19/2021  d_{CDS,PRT}: arithmetic average is used
//                    reinsert: sort by Change::arcDist decsending
// 1.7.1  04/08/2021  inc/indiscern: pairs of indiscernible objects
//                    hybrids are identified for indiscernibility classes
//                    getMissingLeafPairs_ancestors(): only for the member with smallest name of an indiscernibility class
//                    inc/good
// 1.6.2  04/06/2021  makeDistTree: hybrid removal: neighbors are optimized only at the last subgraph optimization iteration
//                                  -delete_deformation_outliers: deformation pair and deformation value are saved
// 1.6.1  03/31/2021  inc/: request2dissim.sh -> pairs2dissim.sh
//                          request_closest.sh -> object2closest.sh
//                          request_closest_sql -> object2closest.sql
//                          grid_min -> pairs2dissim.grid and object2closest.grid
// 1.5.10 03/25/2021  inc/request_closest.sh $OBJ may return $OBJ
// 1.5.9  03/24/2021  NewLeaf::saveRequest() prints empty list of requests if location.indiscernible()
// 1.5.8  03/22/2021  makeDistTree: quality calculation uses threads
// 1.5.7  03/21/2021  FeatureTree::globalSingletonsSize bug
// 1.5.6  03/17/2021  feature tree: only-optional features are removed; 
//                                  Genome::coreNonSingletons does not include optionalCore[]  
//                                  oneFeatureInTree => Genome::coreSet has no optional nominal attributes
//                                  use threads
// 1.5.5  02/28/2021  DTNode::mtx, Subgraph::subPaths2tree() uses threads
// 1.5.4  02/25/2021  combine_dissims.cpp: linear combination
//                    d_symbet
// 1.5.3  02/24/2021  reinsert: arcDist >= max ((size_t) log (name2leaf. size () + 1), areaRadius_std)
// 1.5.2  02/15/2021  PositiveAverageModel::ignoreZero
// 1.5.1  01/31/2021  inc/request_closest.sh has 2 parameters; inc/qc_object.sh
// 1.4.2  01/28/2021  distTree_inc_new_log.sh
// 1.4.1  01/24/2021  inc/{large,request_closest_sql} /dev/null or absent
// 1.3.2  01/10/2021  GeneMark2CDS prints incomplete and ambiguous proteins as well
// 1.3.1  12/18/2020  inc/phen_large --> inc/large; inc/request_closest_sql
// 1.2.3  10/03/2020  makeDistTree does not stop on criterion increase of dissimilarity types > 1
// 1.2.2  08/28/2020  distTree_new.cpp accepts a tree file as input
// 1.2.1  08/20/2020  request2dissim.sh: 4 parameters
// 1.1.12 07/28/2020  Newick import removes transient nodes; printDistTree: 'itree' --> 'dm'
// 1.1.11 07/27/2020  improvement of Func::findZero(); statDistTree: length precision is increased; DistTree::saveFeatureTree(withTime); Feature::getLambda()
// 1.1.10 07/18/2020  makeDistTree -min_arc_prob; DistTree::reroot() does not use 3*SD
// 1.1.9  06/23/2020  DistTree::mergeDissimAttrs() bug affecting indiscernibles
// 1.1.8  06/21/2020  variance_min > 0 => indiscernibles are used, otherwise max_degree is very large
// 1.1.7  06/20/2020  compbine_dissims.cpp: swap(name1,name2)
// 1.1.6  06/17/2020  inc/GenBank/ -> inc/genome/, inc/marker/ -> inc/rRNA/
// 1.1.5  06/16/2020  inc/ contains: server, database, bulk, bulk_remote
// 1.1.4  06/11/2020  distTree_inc_save.sh: inc/*.sh are usable for users
// 1.1.3  06/10/2020  distTree_inc_complete.sh: #3 (update database) parameter is removed
// 1.1.2  06/09/2020  distTree_inc_new.sh: ADD = max($OBJS * $RATE, 1)
// 1.1.1  06/05/2020  [ 1 == 1 ] -> true; [ 1 == 0 ] -> false
// 1.1.0  06/03/2020  inc/request_closest.sh: output is 1 column
// 1.0.2  06/02/2020  cpp/genetics/, cpp/dissim/
// 0.0.0  05/11/2020  align.cpp: DNA: gap_open = -1 (was: 0)
// 0.0.0  05/03/2020  Genome.surveillance = 1: -> dead
// 0.0.0  03/04/2020  uniColl..Genogroup2outliers: "update Genome set genogroup = null" before setting Genome.genogroup
// 0.0.0  02/28/2020  bacteria genomes: BUSCO.10, variance power = 4.0, root at 2780458:5196388 (Thermotogae, not Terabacteria group, cf. 16S tree)
// 0.0.0  02/13/2020  extractFastaProt -replace: print all replacements
// 0.0.0  02/12/2020  hmmsearch2besthits.cpp: print all HMMs (which may hit the same proteins)
// 0.0.0  01/10/2020  update Taxroot set hashes_max = 150000 where id = 4751 and hashes_max = 30000; 
// 0.0.0  11/06/2019  dbo.TaxLineage and phen/ files are added ' [TYPE_SPECIES]' 
// 0.0.0  10/01/2019  distTree_inc.sh: -reinsert is removed
// 0.0.0  09/30/2019  GenBank/: hashed directores: genome, phen; inc/phen_large
// 0.0.0  09/26/2019  panfs/marker/.../inc/request_closest.sh: BLAST output is ordered by nident instead of e-value 
// 0.0.0  09/01/2019  distTree_inc_new.sh: makeDistTree -reinsert 
// 0.0.0  09/01/2019  better optimization for weights \sim 0
// 0.0.0  09/01/2019  Genogroup2outliers: p-value is lowered (two-tail test): Genome.genogroup_outlier is set to null
// 0.0.0  09/01/2019  normalized object criterion; normalized arc error density
// 0.0.0  09/01/2019  hybrid2db.sh, db2unhybrid.sh: obsolete; Genome.outlier = 'auto_hybrid': redo
// 0.0.0  08/20/2019  findHybrids(): only reliable decisions are left; no unhybrid'ing 

