// makeDistTree.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Make a distance tree
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "distTree.hpp"
using namespace DistTree_sp;
#include "../version.inc"



namespace 
{
  
 
const string criterionOutlier_definition ("normalized object criterion");
const string deformationOutlier_definition ("relative object deformation");



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Optimize or modify a least-squares distance tree")
	{
	  version = VERSION;
	  
		// Input
	  addKey ("input_tree", "Directory with a tree of " + dmSuff + "-files ending with '/' or a tree file. If empty then neighbor-joining");

	  addKey ("data", dmSuff + "-file without " + strQuote (dmSuff) + "; or directory with data for an incremental tree ending with '/'");
	  addKey ("dissim_attr", "Dissimilarity attribute name in the <data> file; if all positive two-way attributes must be used then ''");
	  addKey ("weight_attr", "Dissimilarity weight attribute name in the <data> file");

	  addKey ("dissim_coeff", "Coefficient to multiply dissimilarity by (after dissim_power is applied)", "1");
	  addKey ("dissim_power", "Power to raise dissimilarity in", "1");

	  addKey ("variance", "Dissimilarity variance function: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
	  addKey ("variance_power", "Power for -variance pow; >= 0", "NaN");
	  addFlag ("variance_dissim", "Variance is computed off dissimilarities");
	  addFlag ("reinsert_variance_dist", "Variance is computed off tree distances for the reinsert optimization");
	  addKey ("variance_min", "Min. dissimilarity variance", "0");  // ; to be added to the computed variance
	  
	  addKey ("good", "List of objects, which should not be outliers");
	  
	  // Processing
	  addKey ("delete", "Delete leaves whose names are in the indicated file");
	  addFlag ("check_delete", "Check that the names to be deleted actually exist in the tree");  
	  addKey ("keep", "Keep only leaves whose names are in the indicated file by deleting all the other leaves");
	  addFlag ("check_keep", "Check that the names to be kept actually exist in the tree");  

	  addFlag ("optimize", "Optimize topology, arc lengths and re-root");
	//addFlag ("whole", "Optimize whole topology, otherwise by subgraphs of radius " + toString (areaRadius_std));
	  addKey ("subgraph_iter_max", "Max. number of iterations of subgraph optimizations over the whole tree; 0 - unlimited", "0");
	  addFlag ("skip_len", "Skip length-only optimization");
	  addFlag ("reinsert", "Reinsert subtrees before subgraph optimizations; works faster if hybrid objects have been removed");
	//addFlag ("reinsert_orig_weights", "Use original weights in the reinsert optimization");  
	  addFlag ("skip_topology", "Skip topology optimization");	  
	  addFlag ("new_only", "Optimize only new objects in an incremental tree, implies not -optimize");  

	  addFlag ("fix_discernible", "Set the indiscernible flag of objects");
	  addFlag ("fix_transient", "Remove transient nodes (nodes with one child)");
	  	
	  addKey ("delete_criterion_outliers", "Delete outliers by " + criterionOutlier_definition + " and save them in the indicated file");  
	  addKey ("criterion_outlier_num_max", "Max. number of outliers ordered by " + criterionOutlier_definition + " descending to delete; 0 - all", "0");
	  
	  addKey ("delete_deformation_outliers", "Delete outliers by " + deformationOutlier_definition + " and save them in the indicated file");  
	  addKey ("deformation_outlier_num_max", "Max. number of outliers ordered by " + deformationOutlier_definition + " descending to delete; 0 - all", "0");

	  addKey ("hybridness_min", "Min. triangle inequality violation for a hybrid object: d(a,b)/(d(a,x)+d(x,b)), > 1", toString (hybridness_min));
	  addKey ("delete_hybrids", "Find hybrid objects with hybridness > hybridness_min, delete them from the tree and save them in the tab-delimited indicated file. Line format: " + string (PositiveAttr2::hybrid_format));
	//addFlag ("delete_all_hybrids", "Iteratively optimize and delete hybrids until all hybrids are deleted");
	  addKey ("hybrid_parent_pairs", "Save parent pairs of hybrid triangles in the tab-delimited indicated file. Line format: " + string (TriangleParentPair::format));
	  addKey ("dissim_boundary", "Boundary between two merged dissmilarities causing discontinuity", toString (dissim_boundary));

	  addFlag ("reroot", "Re-root");
	  addFlag ("root_topological", "Root minimizes average topologcal depth, otherwise average length to leaves weighted by subtree length");
	  addKey  ("reroot_at", string ("Interior node denoted as \'A") + DistTree::objNameSeparator + "B\', which is the LCA of A and B. Re-root above the LCA in the middle of the arc");

	  addKey ("min_arc_prob", "Min. arc probability to retain", "0");

	  addKey ("dissim_request", "Output file with requests to compute needed dissimilarities, tab-delimited line format: <obj1> <obj2>");

    // Output
	  addFlag ("noqual", "Do not compute quality statistics");
	  addKey ("output_tree", "Save the tesulting tree");
	  addKey ("output_tree_tmp", "Save resulting trees after intermediary steps");
	  addKey ("output_feature_tree", "Resulting tree in feature tree format");
	  addFlag ("feature_tree_time", "Add arc time to <output_feature_tree>");
	  addKey ("output_dissim_coeff", "Save the dissimilarity coefficients for all dissimilarity types");
	  addKey ("leaf_errors", "Output Data Master file without " + strQuote (dmSuff) + " with " + criterionOutlier_definition + " and " + deformationOutlier_definition + " for each leaf");
	  addKey ("arc_existence", "Output Data Master file without " + strQuote (dmSuff) + " with length and existence probability for each arc");
	  addKey ("output_dissim", "Output file with dissimilarities used in the tree, tab-delimited line format: <obj1> <obj2> <dissimilarity>");
	//addFlag ("deredundify_indiscernible", "Remove dissimilarities from the file <output_dissim> for indiscernible objects with non-smallest names in indiscernibility classes");
    addFlag ("output_dist_etc", "Add columns " + string (DistTree::dissimExtra) + " to the file <output_dissim>");
    addKey ("output_data", "Dataset file without " + strQuote (dmSuff) + " with merged dissimilarity attribute and dissimilarity variance");
	}
	
	
	
	void deleteHybrids (DistTree &tree,
	                    bool optimizeLeafNeighbors,
	                    OFStream *triangleParentPairsOs,
	                    OFStream &hybridTrianglesOs,
	                    Vector<Pair<const Leaf*>> *hybridDissimRequests) const
	// Append: *hybridDissimRequests
	{
	#if 1
    tree. setLeafNormCriterion ();
  #else
    // Too few hybrids
    tree. setNodeMaxDeformationDissimNum (); 
  #endif

    const Vector<TriangleParentPair> triangleParentPairs (tree. findHybrids (1e0, hybridDissimRequests));  // PAR 
    if (hybridDissimRequests)
      cout << "# Hybrid dissimilarity requests: " << hybridDissimRequests->size () << endl;

    Vector<Triangle> hybridTriangles;
    for (const TriangleParentPair& tpp : triangleParentPairs)
    {    	
    	if (triangleParentPairsOs)
			  tpp. print (*triangleParentPairsOs);
			tpp. qc ();
		  hybridTriangles << tpp. getHybridTriangles ();  
		}
		if (triangleParentPairsOs)
		  *triangleParentPairsOs << endl;

    hybridTriangles. sort ();
    hybridTriangles. uniq ();

    VectorPtr<Leaf> hybrids;  hybrids. reserve (hybridTriangles. size () * 2);  // PAR
    for (const Triangle& tr : hybridTriangles)
    {
    	tr. print (hybridTrianglesOs);
    	tr. qc ();
			ASSERT (tr. hybridness >= hybridness_min);
			hybrids << tr. getHybrids (true);
    }
    hybrids. sort ();
    hybrids. uniq ();
    
    for (const TriangleParentPair& tpp : triangleParentPairs)
      tpp. qcMatchHybrids (hybrids);
    
    {
      Set<const Leaf*> extra;  // These objects will not be in hybridTriangle's !
      for (const Leaf* leaf : hybrids)
      	if (! leaf->discernible)
      		for (const DiGraph::Node* child : leaf->getParent () -> getChildren ())
      		{
      			const Leaf* other = static_cast <const DTNode*> (child) -> asLeaf ();
      			ASSERT (other);
   			  	ASSERT (! other->discernible);
    			  if (! hybrids. containsFast (other))
    			  	extra << other;
      	  }
      Common_sp::insertAll (hybrids, extra);
    }
    			  	      
    couterr << "# Hybrids: " << hybrids. size () << endl;        

    section ("Deleting hybrids", false);
    {
      Progress prog (hybrids. size ());
      for (const Leaf* leaf : hybrids)
      {
      	ASSERT (leaf->graph);
        tree. removeLeaf (var_cast (leaf), optimizeLeafNeighbors);
        prog (tree. absCriterion2str ());
      }
    }
    
    tree. reportErrors (cout);
  //cout << endl;
    tree. qc ();
	}
	
	
	
	void body () const final
  {
	  const string input_tree          = getArg ("input_tree");
	  const string dataFName           = getArg ("data");
	  const string dissimAttrName      = getArg ("dissim_attr");
	  const string multAttrName        = getArg ("weight_attr");
	               dissim_power        = str2real (getArg ("dissim_power"));      // Global
	               dissim_coeff        = str2real (getArg ("dissim_coeff"));      // Global
	               varianceType        = str2varianceType (getArg ("variance"));  // Global
	               variancePower       = str2real (getArg ("variance_power"));    // Global
	               variance_min        = str2real (getArg ("variance_min"));      // Global
	  const bool   variance_dissim     = getFlag ("variance_dissim");
	  const bool   reinsert_variance_dist = getFlag ("reinsert_variance_dist");
	  const string goodFName           = getArg ("good");
	               
		const string deleteFName         = getArg ("delete");
		const bool   check_delete        = getFlag ("check_delete");
		const string keepFName           = getArg ("keep");
		const bool   check_keep          = getFlag ("check_keep");
		      
		const bool   optimize            = getFlag ("optimize");
	//const bool   whole               = getFlag ("whole");
		const size_t subgraph_iter_max   = str2<size_t> (getArg ("subgraph_iter_max"));
		const bool   skip_len            = getFlag ("skip_len");
		const bool   reinsert            = getFlag ("reinsert");
	//const bool   reinsert_orig_weights = getFlag ("reinsert_orig_weights");		
		const bool   skip_topology       = getFlag ("skip_topology");
	  const bool   new_only            = getFlag ("new_only");

	  const bool   fix_discernible     = getFlag ("fix_discernible");
	  const bool   fix_transient       = getFlag ("fix_transient");
		
		const bool   reroot              = getFlag ("reroot");		
		const bool   root_topological    = getFlag ("root_topological");
		const string reroot_at           = getArg ("reroot_at");

  	const Prob   arcExistence_min    = str2real (getArg ("min_arc_prob"));

		const string delete_criterion_outliers = getArg ("delete_criterion_outliers");
		const size_t criterion_outlier_num_max = str2<size_t> (getArg ("criterion_outlier_num_max"));

		const string delete_deformation_outliers = getArg ("delete_deformation_outliers");
		const size_t deformation_outlier_num_max = str2<size_t> (getArg ("deformation_outlier_num_max"));

		const string delete_hybrids      = getArg ("delete_hybrids");
	//const bool   delete_all_hybrids  = getFlag ("delete_all_hybrids");
		const string hybrid_parent_pairs = getArg ("hybrid_parent_pairs");
                 dissim_boundary     = str2real (getArg ("dissim_boundary"));
		             hybridness_min      = str2real (getArg ("hybridness_min"));
		
		const bool   noqual              = getFlag ("noqual");

		const string output_tree         = getArg ("output_tree");
		const string output_tree_tmp     = getArg ("output_tree_tmp");
		const string output_dissim_coeff = getArg ("output_dissim_coeff");
		const string output_feature_tree = getArg ("output_feature_tree");
		const bool   feature_tree_time   = getFlag ("feature_tree_time");
		const string leaf_errors         = getArg ("leaf_errors");
		const string arc_existence       = getArg ("arc_existence");
		const string dissim_request      = getArg ("dissim_request");
		const string output_dissim       = getArg ("output_dissim");
  //const bool   deredundify_indiscernible = getFlag ("deredundify_indiscernible");
		const bool   output_dist_etc     = getFlag ("output_dist_etc");
		const string output_data         = getArg ("output_data");
		
		const bool optimizable = isDirName (dataFName) || ! dataFName. empty ();
		
    if (isDirName (dataFName))
    {
      if (! dissimAttrName. empty ())
        throw runtime_error ("Non-empty dissimilarity attribute with no " + dmSuff + "-file");
      if (! multAttrName. empty ())
        throw runtime_error ("Non-empty dissimilarity weight attribute with no " + dmSuff + "-file");
    }
    else
      if (dataFName. empty () && ! dissimAttrName. empty ())
        throw runtime_error ("Dissimilarity attribute with no data file");
    if (dissimAttrName. empty () && ! multAttrName. empty ())
      throw runtime_error ("Dissimilarity weight attribute with no dissimilarity attribute");

    if (dissim_coeff <= 0.0)
      throw runtime_error ("-dissim_coeff must be positive");
    if (dissim_power <= 0.0)
      throw runtime_error ("-dissim_power must be positive");

    if (variance_min < 0.0)
      throw runtime_error ("-variance_min cannot be negative");
    if (variance_min && varianceType == varianceType_none)
      throw runtime_error ("-variance_min requires a variance function");
    if (varianceType != varianceType_none && dataFName. empty ())
      throw runtime_error ("Variance function requires a data file");
    if (varianceType != varianceType_none && ! multAttrName. empty ())
      throw runtime_error ("Variance function excludes a dissimilarity weight attribute");
		if (! isNan (variancePower) && varianceType != varianceType_pow) 
		  throw runtime_error ("-variance_power requires -variance pow");
		if (isNan (variancePower) && varianceType == varianceType_pow)
		  throw runtime_error ("-variance pow requires -variance_power");
		if (variancePower < 0.0)
		  throw runtime_error ("-variance_power must be non-negative");

    if (check_delete && deleteFName. empty ())
      throw runtime_error ("-check_delete requires -delete");
    if (check_keep && keepFName. empty ())
      throw runtime_error ("-check_keep requires -keep");
    if (! deleteFName. empty () && ! keepFName. empty ())
      throw runtime_error ("Cannot use both -delete and -keep");

    if (optimize && ! optimizable)
      throw runtime_error ("-optimize requires dissimilarities");
    if (subgraph_iter_max && ! optimize)
      throw runtime_error ("-subgraph_iter_max requires -optimize");
    if (skip_len && ! optimize)
      throw runtime_error ("-skip_len requires -optimize");
    if (reinsert && ! optimize)
      throw runtime_error ("-reinsert requires -optimize");
  //if (reinsert_orig_weights && ! reinsert)
    //throw runtime_error ("-reinsert_orig_weights requires -reinsert");
    if (skip_topology && ! optimize)
      throw runtime_error ("-skip_topology requires -optimize");
    if (new_only && optimize)
      throw runtime_error ("-new_only excludes -optimize");
      
    if (fix_discernible && dataFName. empty ())
      throw runtime_error ("-fix_discernible requires dissimilarities");
    if (fix_discernible && variance_min)
      throw runtime_error ("-fix_discernible requires zero -variance_min");
      
    if (! delete_hybrids. empty () && ! optimizable)
      throw runtime_error ("-delete_hybrids requires dissimilarities");
    if (dissim_boundary <= 0.0)
    	throw runtime_error ("-dissim_boundary must be > 0");
    if (hybridness_min <= 1.0)
    	throw runtime_error ("-hybridness_min must be > 1");
  //if (delete_all_hybrids && delete_hybrids. empty ())
    //throw runtime_error ("-delete_all_hybrids requires -delete_hybrids");
    if (! hybrid_parent_pairs. empty () && delete_hybrids. empty ())      
    	throw runtime_error ("-hybrid_parent_pairs requires -delete_hybrids");

    if (! delete_criterion_outliers. empty () && ! optimizable)
      throw runtime_error ("-delete_criterion_outliers requires dissimilarities");
    if (criterion_outlier_num_max && delete_criterion_outliers. empty ())
      throw runtime_error ("-criterion_outlier_num_max requires -delete_criterion_outliers");
    if (deformation_outlier_num_max && delete_deformation_outliers. empty ())
      throw runtime_error ("-deformation_outlier_num_max requires -delete_deformation_outliers");

		if (reroot && ! reroot_at. empty ())
		  throw runtime_error ("-reroot excludes -reroot_at");

		if (arcExistence_min && ! optimizable)
		  throw runtime_error ("-min_arc_prob requires dissimilarities");
		if (! isProb (arcExistence_min))
		  throw runtime_error ("-min_arc_prob must be between 0 and 1");
		  
    if (! output_dissim. empty () && ! optimizable)
      throw runtime_error ("-output_dissim requires dissimilarities");    	
    if (output_dist_etc && output_dissim. empty ())
    	throw runtime_error ("-output_dist_etc requires -output_dissim");
  //if (deredundify_indiscernible && output_dissim. empty ())
    //throw runtime_error ("-deredundify_indiscernible requires -output_dissim");
    if (! leaf_errors. empty () && ! optimizable)
      throw runtime_error ("-leaf_errors requires dissimilarities");    	
    if (! leaf_errors. empty () && noqual)
      throw runtime_error ("-noqual excludes -leaf_errors");
    if (! dissim_request. empty () && ! optimizable)
      throw runtime_error ("-dissim_request requires dissimilarities");    	
    if (! output_data. empty () && ! optimizable)
      throw runtime_error ("-output_data requires dissimilarities");    	



    DistTree::printParam (cout);
    cout << "Root: " << (root_topological ? "topological" : "by length") << endl;
    cout << endl;


    unique_ptr<DistTree> tree;
    {
      const Chronometer_OnePass cop ("Initial topology");  
      tree. reset (isDirName (dataFName)
                     ? new DistTree (dataFName, input_tree, true, true, new_only)
                     : input_tree. empty ()
                       ? new DistTree (            dataFName, dissimAttrName, multAttrName)
                       : new DistTree (input_tree, dataFName, dissimAttrName, multAttrName)
                  );
    }
    ASSERT (tree. get ());
    if (tree->getDiscernibles (). size () == 1)
    {
      cout << "One discernible object" << endl;
      return;
    }
    
    ASSERT (optimizable == tree->optimizable ());
    
    if (variance_dissim)
    {
      if (tree->multFixed)
        throw runtime_error ("-variance_dissim cannot be used with fixed dissimilarity variance");
      else
      {
      	tree->setDissimMult (false);  
        tree->multFixed = true;
      }
    }
    
    if (fix_discernible)
    {
      section ("Fixing discernible", false);
      tree->setDiscernibles ();
      const Keep<bool> multFixed_old (tree->multFixed);
      tree->multFixed = false;
      tree->setDissimMult (! multFixed_old. get ());
    }
    
    if (fix_transient)
      tree->fixTransients ();
      
    if (! goodFName. empty ())
      tree->setGoodLeaves (goodFName);
      
    if (qc_on)
    {
      section ("QC", false);
      tree->qc (); 
    }

    tree->printInput (cout);
    cout << endl;
    
    
    
    if (! deleteFName. empty ())
    {
      section ("Deleting", false);
      size_t deleted = 0;
      {
        LineInput f (deleteFName);  
        Progress prog (0, 1000);  // PAR
        while (f. nextLine ())
        {
          trim (f. line);
          const Leaf* leaf = findPtr (tree->name2leaf, f. line);
          if (! leaf)
          {
            if (check_delete)
              throw runtime_error ("Leaf " + f. line + " not found");
            continue;
          }
          tree->removeLeaf (var_cast (leaf), false);
         	prog ();
          deleted++;
        }
      }
      cout << "# Deleted: " << deleted << endl;
	    if (tree->optimizable ())
        tree->reportErrors (cout);
      tree->qc ();
      cout << endl;
    }


    if (! keepFName. empty ())
    {
      section ("Keeping", false);
      VectorPtr<Leaf> toDelete;  toDelete. reserve (tree->name2leaf. size ());
      for (const auto& it : tree->name2leaf)
        toDelete << it. second;
      toDelete. sort ();
      ASSERT (toDelete. isUniq ());
      VectorPtr<Leaf> toKeep;  toKeep. reserve (tree->name2leaf. size ());
      {
        LineInput f (keepFName, 1000);  // PAR
        while (f. nextLine ())
        {
          trim (f. line);
          const Leaf* leaf = findPtr (tree->name2leaf, f. line);
          if (! leaf)
          {
            if (check_keep)
              throw runtime_error ("Leaf " + f. line + " not found");
            continue;
          }
          toKeep << leaf;  
        }
      }
      toKeep. sort ();
      toKeep. uniq ();
      toDelete. setMinus (toKeep);
      {
        Progress prog (0, 1000);  // PAR
        for (const Leaf* leaf : toDelete)
        {
          tree->removeLeaf (var_cast (leaf), false);
          if (tree->optimizable ())
            prog (tree->absCriterion2str ());
          else
          	prog ();
        }
      }
      cout << "# Deleted: " << toDelete. size () << endl;
	    if (tree->optimizable ())
        tree->reportErrors (cout);
      tree->qc ();
      cout << endl;
    }


    Vector<Pair<const Leaf*>> hybridDissimRequests;
    unique_ptr<OFStream> hybridParentPairsF;
    if (! hybrid_parent_pairs. empty ())
    	hybridParentPairsF. reset (new OFStream (hybrid_parent_pairs));
    unique_ptr<OFStream> hybridF;
    if (! delete_hybrids. empty ())
    	hybridF. reset (new OFStream (delete_hybrids));
    if (tree->optimizable ())
    {
      if (optimize)
      {
        const size_t leaves = tree->root->getLeavesSize ();
        if (leaves > 3)
        {
          tree->saveFile (output_tree_tmp);  

          bool predictionImproved = false;

          if (! skip_len)            
          {
            const Chronometer_OnePass cop ("Initial arc lengths");

            section ("Optimizing topology: arc lengths for the whole tree", true);
            tree->optimizeLenWhole ();
            cout << tree->absCriterion2str () << endl;
            tree->saveFile (output_tree_tmp); 

            section ("Optimizing topology: arc lengths at each arc", true);
            const size_t lenArc_deleted = tree->optimizeLenArc ();
            cout << "# Nodes deleted = " << lenArc_deleted << endl;
            cout << tree->absCriterion2str () << endl;
            tree->saveFile (output_tree_tmp); 

            section ("Optimizing topology: arc lengths at each node", true);
            const size_t lenNode_deleted = tree->optimizeLenNode ();
            cout << "# Nodes deleted = " << lenNode_deleted << endl;
            cout << tree->absCriterion2str () << endl;
            tree->saveFile (output_tree_tmp); 
            
            tree->qc ();
            if (verbose ())
              tree->saveText (cout);  
            tree->reportErrors (cout);
            
            if (lenArc_deleted || lenNode_deleted)
              predictionImproved = true;
          }
          

          if (reinsert)
          {
            section ("Optimizing topology: reinsert", true);
            const Chronometer_OnePass cop ("Topology optimization: reinsert");
            if (! tree->multFixed && ! reinsert_variance_dist)
            {
         	    tree->setDissimMult (false);   // may damage topology for big variance functions
              couterr << tree->absCriterion2str () << endl; 
         	  }
            tree->optimizeReinsert ();  
            tree->saveFile (output_tree_tmp); 
            predictionImproved = true;
          }
          
          if (predictionImproved)
          {
         		tree->setDissimMult (true);
            couterr << tree->absCriterion2str () << endl; 
         	}
          
          if (! skip_topology)
          {
            section ("Optimizing topology: subgraphs", true);
            const Chronometer_OnePass cop ("Topology optimization: local");
            size_t iter_max = numeric_limits<size_t>::max ();
            if (subgraph_iter_max)
            	minimize (iter_max, subgraph_iter_max);
            ASSERT (iter_max);
            size_t iter = 0;
            while (iter < iter_max)
          	{
              section ("Iteration " + to_string (iter + 1) + ifS (iter_max < numeric_limits<size_t>::max (), " / " + to_string (iter_max)), true);
          		const Real absCriterion_old = tree->absCriterion;
          		ASSERT (absCriterion_old < inf);
              tree->optimizeDissimCoeffs ();  
              tree->optimizeLargeSubgraphs (nullptr);
              cout << tree->absCriterion2str () << endl; 
              if (hybridF. get ())
              	deleteHybrids ( *tree
              	              , iter + 1 == iter_max
              	              , hybridParentPairsF. get ()
              	              , *hybridF
              	              , dissim_request. empty () ? nullptr : & hybridDissimRequests
              	              );
              tree->saveFile (output_tree_tmp); 
              iter++;
              if (! absCriterion_old)
              	break;
              if (   absCriterion_old <= tree->absCriterion 
                  && tree->dissimTypesNum () == 1
                 )
              	break;
              tree->setDissimMult (true);
              if (! tree->multFixed)
                couterr << tree->absCriterion2str () << endl; 
            }
            cout << "# Iterations of subgraph optimization: " << iter << endl;
            tree->reportErrors (cout);
          }
          
          section ("Re-rooting", false);
          const Real radius_ave = tree->reroot (root_topological);
				  const ONumber on (cout, dissimDecimals / 2, false);  // PAR
          cout << "Ave. radius: " << radius_ave << endl;
          cout << endl;
        }
        else if (leaves == 3)
          tree->optimize3 ();
        else if (leaves == 2)
          tree->optimize2 ();
        else
          { ERROR; }
      }
      else
        if (hybridF. get ())
        {
          deleteHybrids (*tree, true, hybridParentPairsF. get (), *hybridF, dissim_request. empty () ? nullptr : & hybridDissimRequests);
          cout << endl;
        }


      if (! delete_criterion_outliers. empty ())
      {
        section ("Finding criterion outliers", false);
	      tree->setLeafNormCriterion (); 
        const Dataset leafErrorDs (tree->getLeafErrorDataset (true, NaN));
	      Real outlier_min_excl = NaN;
	      const VectorPtr<Leaf> outliers (tree->findCriterionOutliers (leafErrorDs, 1e-6, outlier_min_excl));  // PAR
	      const ONumber on (cout, absCriterionDecimals, false);  
	      cout << "# Criterion outliers: " << outliers. size () << endl;
	      cout << "# Criterion outlier threshold: " << outlier_min_excl << endl;
        OFStream f (delete_criterion_outliers);
	      if (! outliers. empty ())
	      {
          section ("Deleting criterion outliers", false);
          {
            Progress prog (outliers. size ());
            size_t removed = 0;
            for (const Leaf* leaf : outliers)
            {
              if (criterion_outlier_num_max && removed >= criterion_outlier_num_max)
                break;
              f << leaf->name << endl;
              tree->removeLeaf (var_cast (leaf), true);
              prog (tree->absCriterion2str ());
              removed++;
            }
            if (removed)
              tree->setDissimMult (true);
          }
          tree->reportErrors (cout);
        }
        cout << endl;
      }
      tree->qc ();
    
      
      if (! delete_deformation_outliers. empty ())
      {
        section ("Finding deformation outliers", false);
	      tree->setNodeMaxDeformationDissimNum (); 
	      Real outlier_min_excl = NaN;
	      const VectorPtr<Leaf> outliers (tree->findDeformationOutliers (tree->getDeformation_mean (), 1e-10, outlier_min_excl));  // PAR  
	      const ONumber on (cout, absCriterionDecimals, false);  
	      cout << "# Deformation outliers: " << outliers. size () << endl;
	      cout << "# Deformation outlier threshold: " << outlier_min_excl << endl;
        OFStream f (delete_deformation_outliers);
	      if (! outliers. empty ())
	      {
          section ("Deleting deformation outliers", false);
          {
            Progress prog (outliers. size ());
            size_t removed = 0;
            for (const Leaf* leaf : outliers)
            {
              if (deformation_outlier_num_max && removed >= deformation_outlier_num_max)
                break;
              f << leaf->name << '\t' << leaf->getDeformationS () << endl;
              tree->removeLeaf (var_cast (leaf), true);
              prog (tree->absCriterion2str ());
              removed++;
            }
            if (removed)
              tree->setDissimMult (true);
          }
          tree->reportErrors (cout);
        }
        for (DiGraph::Node* node : tree->nodes)
          static_cast <DTNode*> (node) -> maxDeformationDissimNum = dissims_max;
        cout << endl;
      }
      tree->qc ();


      if (! noqual)
      {
        section ("Node/arc criteria", true);
        tree->setLeafNormCriterion ();  
        tree->setNodeMaxDeformationDissimNum ();
        tree->setErrorDensities ();  
        tree->qc ();
      }


      cout << "OUTPUT:" << endl;  
      tree->reportErrors (cout);
      cout << endl;
    }
    
    
    chron_getBestChange.    print (cout);
    chron_tree2subgraph.    print (cout);
    chron_subgraphOptimize. print (cout); 
    chron_subgraph2tree.    print (cout); 


    // Tree model is fixed
    

    if (reroot)
    {
      section ("Re-rooting", false);
      cout << "Ave. radius: " << tree->reroot (root_topological) << endl;
      cout << endl;
      tree->qc ();
    }
    else 
      if (! reroot_at. empty ())
      {
			  Tree::LcaBuffer buf;
        const DTNode* underRoot = tree->lcaName2node (reroot_at, buf);
        tree->reroot (var_cast (underRoot), underRoot->len / 2.0);
        tree->qc ();
      }
      

    if (arcExistence_min)
    {
      const size_t arcsDeleted = tree->deleteQuestionableArcs (arcExistence_min);
      cout << "# Arcs deleted: " << arcsDeleted << endl;
      cout << endl;
    }


    if (! noqual && tree->optimizable ())
    {
		  const ONumber on (cout, 2, false);  // PAR
    	const Real dissim_var = tree->getUnoptimizable ();  
      cout << "Relative epsilon2_0 = " << sqrt (dissim_var / tree->target2_sum) * 100.0 << " %" << endl;
        // Must be << "Average arc error"
      cout << "Mean residual = " << tree->getMeanResidual () << endl;
      cout << "Correlation between residual^2 and dissimilarity = " << tree->getSqrResidualCorr () << endl; 
      cout << endl;
      tree->qc ();
    }
    
    
    tree->saveFile (output_tree); 
    tree->saveDissimCoeffs (output_dissim_coeff);
    tree->saveFeatureTree (output_feature_tree, feature_tree_time);


    if (! output_data. empty ())
    {
      Real dissimTypeError = NaN;
      {
        tree->dissims. sort ();
        const Dataset ds (tree->getDissimWeightDataset (dissimTypeError));
        OFStream of (output_data + dmSuff);
        ds. saveText (of);
      }
      if (! tree->dissimTypes. empty ())
      {
        const ONumber on (cout, absCriterionDecimals, false);  
        cout << "Error between dissimilarities of different types = " << dissimTypeError << endl;
        tree->reportErrors (cout, dissimTypeError);
        cout << endl;
      }
    }
    
    
    if (! leaf_errors. empty ())
    {
      tree->setLeafNormCriterion (); 
      tree->setNodeMaxDeformationDissimNum ();
      const Dataset leafErrorDs (tree->getLeafErrorDataset (true, tree->getDeformation_mean ()));
      OFStream f (leaf_errors + dmSuff);
	    leafErrorDs. saveText (f);    
    }


    if (! arc_existence. empty ())
    {
      section ("Arc existence", false);
      Dataset ds;
      auto len = new PositiveAttr1 ("len", ds, dissimDecimals);
      auto prob = new ProbAttr1 ("prob", ds, 3);  // PAR
      Real realArcs = 0.0;
      size_t n = 0;
      {
        const size_t nodesSize = tree->nodes. size ();
        Progress prog (nodesSize, 1000);  // PAR
        for (const DiGraph::Node* node : tree->nodes)
        {
          prog ();
          const DTNode* dtNode = static_cast <const DTNode*> (node);
          if (dtNode == tree->root)
            continue;
          if (dtNode->inDiscernible ())
            continue;
          const size_t objNum = ds. appendObj (dtNode->getLcaName ());
          const Prob arcExists = dtNode->getArcExistence ();
          (*len)  [objNum] = dtNode->len;
          (*prob) [objNum] = arcExists;
          realArcs += arcExists;
          n++;
        }
      }
      ds. qc ();
      OFStream of (arc_existence + dmSuff);
      ds. saveText (of);
      const ONumber on (cout, 3, false);  // PAR
      cout << "Fraction of real interior arcs: " << realArcs / (Real) n << " (" << realArcs << " / " << n << ')' << endl;
    }


    if (! output_dissim. empty ())
    {
      OFStream f (output_dissim);
      tree->saveDissim (f, true /*! deredundify_indiscernible*/, output_dist_etc);
    }
    
    
    if (! dissim_request. empty ())
    {
      section ("Finding missing leaf pairs", false);
      Vector<Pair<const Leaf*>> pairs (tree->getMissingLeafPairs_ancestors (sparsingDepth, false));
      cout << endl;
      cout << "# Ancestor-based dissimilarity requests: " << pairs. size () << endl;

    #if 0
      // Needed if closest objects are not available
      ??
      {      
        const VectorPtr<Leaf> depthOutliers (tree->findDepthOutliers ());
        const Vector<Pair<const Leaf*>> depthPairs (tree->leaves2missingLeafPairs (depthOutliers));
        cout << "# Depth-outlier-based requests: " << depthPairs. size () << endl;
        pairs << depthPairs;  // --> unionFast () ??
      }
      
      {
        // Clusters of diameter <= diameter_max ??
        const VectorPtr<DTNode> depthClusters (tree->findDepthClusters (100));  // PAR
      //cout << "# Depth clusters: " << depthClusters. size () << endl;
        VectorPtr<Leaf> leaves;  leaves. reserve (depthClusters. size () * 4);
        for (const DTNode* node : depthClusters)
        {
          leaves << node->reprLeaf;
        	for (const DiGraph::Arc* arc : node->arcs [false])
        	  leaves <<  static_cast <const DTNode*> (arc->node [false]) -> reprLeaf;
        }
        leaves. sort ();
        leaves. uniq ();
        cout << "# Depth leaves: " << leaves. size () << endl;
        const Vector<Pair<const Leaf*>> depthPairs (tree->leaves2missingLeafPairs (leaves));
        cout << "# Depth-based requests: " << depthPairs. size () << endl;
        pairs << depthPairs;  // --> unionFast () ??        
      }
    #endif
    
      pairs << hybridDissimRequests;      
      pairs. sort ();
      pairs. uniq ();            
      cout << "# Total dissimilarity requests: " << pairs. size () << endl;
      {
        OFStream f (dissim_request);      
        for (const auto& p : pairs)
          f << p. first->name << '\t' << p. second->name << endl;
      }
    }
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


