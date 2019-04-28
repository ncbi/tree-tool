// makeDistTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{
  


void checkOptimizable (const DistTree& tree,
                       const string &parameter)
{
  if (! tree. optimizable ())
    throw runtime_error ("Parameter " + parameter + " requires dissimilarities");
}



const string outlierCriterion ("relative unweighted average leaf error");



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Optimize topology and compute statistics for a least-squares distance tree")
	{
		// Input
	  addKey ("input_tree", "Directory with a tree of " + dmSuff + "-files ending with '/' or a tree file. If empty then neighbor-joining");
	  addKey ("data", dmSuff + "-file without " + strQuote (dmSuff) + "; or directory with data for an incremental tree ending with '/'");
	  addKey ("dissim_attr", "Dissimilarity attribute name in the <data> file; if all positive two-way attributes must be used then ''");
	  addKey ("weight_attr", "Dissimilarity weight attribute name in the <data> file");
	  addKey ("dissim_power", "Power to raise dissimilarity in", "1");
	  addKey ("dissim_coeff", "Coefficient to multiply dissimilarity by (after dissim_power is applied)", "1");
	  addKey ("variance", "Dissimilarity variance function: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
	  addKey ("var_min", "Min. dissimilarity variance; to be added to the computed variance. If > 0 then all objects are discernible", "0");
	  addFlag ("refresh_dissim", "Add more requests to <dissim_request>");
	  
	  // Processing
	  addKey ("delete", "Delete leaves whose names are in the indicated file");
	  addFlag ("check_delete", "Check that the names to be deleted actually exist in the tree");  
	  addKey ("keep", "Keep only leaves whose names are in the indicated file by deletign all the other leaves");
	  addFlag ("check_keep", "Check that the names to be kept actually exist in the tree");  

	  addFlag ("optimize", "Optimize topology, arc lengths and re-root");
	//addFlag ("whole", "Optimize whole topology, otherwise by subgraphs of radius " + toString (areaRadius_std));
	  addFlag ("subgraph_fast", "Limit the number of iterations of subgraph optimizations over the whole tree");
	  addKey ("subgraph_iter_max", "Max. number of iterations of subgraph optimizations over the whole tree; 0 - unlimited", "0");
	  addFlag ("skip_len", "Skip length-only optimization");
	  addFlag ("reinsert", "Re-insert subtrees before subgraph optimizations");
	  addFlag ("reinsert_iter", "Re-insert subtrees at each iteration of subgraph optimization");
	  addFlag ("skip_topology", "Skip topology optimization");	  
	  addFlag ("new_only", "Optimize only new objects in an incremental tree, implies not -optimize");  
	  addFlag ("fix_discernibles", "Set the indiscernible flag of objects");

	  addFlag ("reroot", "Re-root");
	  addFlag ("root_topological", "Root minimizes average topologcal depth, otherwise average length to leaves weighted by subtree length");
	  addKey  ("reroot_at", string ("Interior node denoted as \'A") + DistTree::objNameSeparator + "B\', which is the LCA of A and B. Re-root above the LCA in the middle of the arc");
	  	
	  addKey ("delete_outliers", "Delete outliers by " + outlierCriterion + " and save them in the indicated file");  
	  addKey ("outlier_num_max", "Max. number of outliers ordered by " + outlierCriterion + " descending to delete; 0 - all", "0");
	  
	  addKey ("delete_hybrids", "Find hybrid objects with hybridness > hybridness_min, delete them from the tree and save them in the tab-delimited indicated file. Line format: " + string (Triangle::format));
	//addFlag ("delete_all_hybrids", "Iteratively optimize and delete hybrids until all hybrids are deleted");
	  addKey ("hybrid_parent_pairs", "Save parent pairs of hybrid triangles in the tab-delimited indicated file. Line format: " + string (TriangleParentPair::format));
	  addKey ("dissim_boundary", "Boundary between two merged dissmilarity measures causing discontinuity", toString (dissim_boundary));
	  addKey ("hybridness_min", "Min. triangle inequality violation for a hybrid object: d(a,b)/(d(a,x)+d(x,b)), > 1", toString (hybridness_min));

	  addFlag ("noqual", "Do not compute quality statistics");

    // Output
	  addKey ("output_tree", "Resulting tree");
	  addKey ("output_dissim_coeff", "Save the dissimilarity coefficients for all dissimilarity types");
	  addKey ("output_feature_tree", "Resulting tree in feature tree format");
	  addKey ("leaf_errors", "Output " + dmSuff + "-file without " + strQuote (dmSuff) + " with " + outlierCriterion + " for each leaf");
	  addKey ("dissim_request", "Output file with requests to compute needed dissimilarities, tab-delimited line format: <obj1> <obj2>");
	  addKey ("output_dissim", "Output file with dissimilarities used in the tree, tab-delimited line format: <obj1> <obj2> <dissim>");
    addFlag ("output_dist_etc", "Add columns " + string (DistTree::dissimExtra) + " to the file <output_dissim>");
    addKey ("output_data", "Dataset file without " + strQuote (dmSuff) + " with merged dissimilarity attribute and dissimilarity variance");
	}
	
	
	
	bool deleteHybrids (DistTree &tree,
	                    bool searchBadLeaves,
	                    OFStream *triangleParentPairsOs,
	                    OFStream &hybridTrianglesOs,
	                    Vector<Pair<const Leaf*>> *hybridDissimRequests) const
	// Return: true <=> hybrids are deleted
	// Append: *hybridDissimRequests
	{
	  if (searchBadLeaves)
      tree. setLeafAbsCriterion ();

    const Vector<TriangleParentPair> triangleParentPairs (tree. findHybrids (1e0, searchBadLeaves, hybridDissimRequests));  // PAR 
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
    
  #if 0
    // Should be in hybridTrianglesOs
    Set<const Leaf*> extra;
    for (const Leaf* leaf : hybrids)
    	if (! leaf->discernible)
    		for (const DiGraph::Node* child : leaf->getParent () -> getChildren ())
    			if (const Leaf* other = static_cast <const DTNode*> (child) -> asLeaf ())
    			  if (! hybrids. containsFast (other))
    			  {
    			  	ASSERT (! other->discernible);
    			  	extra << other;
    			  }
    Common_sp::insertAll (hybrids, extra);
  #endif
    			  	      
    cout << "# Hybrids: " << hybrids. size () << endl;        

    cerr << "Deleting hybrids ..." << endl;
    bool deleted = false;
    {
      Progress prog (hybrids. size ());
      for (const Leaf* leaf : hybrids)
      {
      	ASSERT (leaf->graph);
        tree. removeLeaf (var_cast (leaf), true);
        prog (tree. absCriterion2str ());
        deleted = true;
      }
    }
    
    tree. reportErrors (cout);
    cout << endl;
    tree. qc ();
    
    return deleted;
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
	               variance_min        = str2real (getArg ("var_min"));
	               
		const string deleteFName         = getArg ("delete");
		const bool   check_delete        = getFlag ("check_delete");
		const string keepFName           = getArg ("keep");
		const bool   check_keep          = getFlag ("check_keep");
		      
		const bool   optimize            = getFlag ("optimize");
	//const bool   whole               = getFlag ("whole");
		const bool   subgraph_fast       = getFlag ("subgraph_fast");
		const size_t subgraph_iter_max   = str2<size_t> (getArg ("subgraph_iter_max"));
		const bool   skip_len            = getFlag ("skip_len");
		const bool   reinsert            = getFlag ("reinsert");
		const bool   reinsert_iter       = getFlag ("reinsert_iter");
		const bool   skip_topology       = getFlag ("skip_topology");
	  const bool   new_only            = getFlag ("new_only");
	  const bool   fix_discernibles    = getFlag ("fix_discernibles");
		
		const bool   reroot              = getFlag ("reroot");		
		const bool   root_topological    = getFlag ("root_topological");
		const string reroot_at           = getArg ("reroot_at");

		const string delete_outliers     = getArg ("delete_outliers");
		const size_t outlier_num_max     = str2<size_t> (getArg ("outlier_num_max"));

                 dissim_boundary     = str2real (getArg ("dissim_boundary"));
		             hybridness_min      = str2real (getArg ("hybridness_min"));
		const string delete_hybrids      = getArg ("delete_hybrids");
	//const bool   delete_all_hybrids  = getFlag ("delete_all_hybrids");
		const string hybrid_parent_pairs = getArg ("hybrid_parent_pairs");
		
		const bool   noqual              = getFlag ("noqual");

		const string output_tree         = getArg ("output_tree");
		const string output_dissim_coeff = getArg ("output_dissim_coeff");
		const string output_feature_tree = getArg ("output_feature_tree");
		const string leaf_errors         = getArg ("leaf_errors");
		const string dissim_request      = getArg ("dissim_request");
		const bool   refresh_dissim      = getFlag ("refresh_dissim");
		const string output_dissim       = getArg ("output_dissim");
		const bool   output_dist_etc     = getFlag ("output_dist_etc");
		const string output_data         = getArg ("output_data");
		
		ASSERT (! (reroot && ! reroot_at. empty ()));
    if (isDirName (dataFName))
    {
      if (! input_tree. empty ())
        throw runtime_error ("Input tree must be in " + dataFName);
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
      throw runtime_error ("dissim_coeff must be positive");
    if (dissim_power <= 0.0)
      throw runtime_error ("dissim_power must be positive");
    if (check_delete && deleteFName. empty ())
      throw runtime_error ("-check_delete requires -delete");
    if (check_keep && keepFName. empty ())
      throw runtime_error ("-check_keep requires -keep");
    if (! deleteFName. empty () && ! keepFName. empty ())
      throw runtime_error ("Cannot use both -delete and -keep");
    if (variance_min < 0.0)
      throw runtime_error ("-var_min cannot be negative");
    if (variance_min && varianceType == varianceType_none)
      throw runtime_error ("-var_min needs a variance function");
    if (varianceType != varianceType_none && dataFName. empty ())
      throw runtime_error ("Variance function needs a data file");
    if (varianceType != varianceType_none && ! multAttrName. empty ())
      throw runtime_error ("Variance function cannot be used with a dissimilarity weight attribute");
    if (output_dist_etc && output_dissim. empty ())
    	throw runtime_error ("output_dist_etc exists, but no output_dissim");
    IMPLY (subgraph_fast,     optimize);
    if (subgraph_iter_max && ! optimize)
      throw runtime_error ("-subgraph_iter_max requires -optimize");
    IMPLY (skip_len,          optimize);
    IMPLY (reinsert,          optimize);
    IMPLY (reinsert_iter,     optimize);
    IMPLY (skip_topology,     optimize);
    if (new_only && optimize)
      throw runtime_error ("-new_only cannot be used with -optimize");
    if (fix_discernibles && dataFName. empty ())
      throw runtime_error ("-fix_discernibles needs dissimilarities");
    if (fix_discernibles && variance_min)
      throw runtime_error ("-fix_discernibles implies zero -var_min");
    if (dissim_boundary <= 0.0)
    	throw runtime_error ("dissim_boundary must be > 0");
    if (hybridness_min <= 1.0)
    	throw runtime_error ("hybridness_min must be > 1");
  //if (delete_all_hybrids && delete_hybrids. empty ())
    //throw runtime_error ("-delete_all_hybrids assumes -delete_hybrids");
    if (! hybrid_parent_pairs. empty () && delete_hybrids. empty ())      
    	throw runtime_error ("-hybrid_parent_pairs assumes -delete_hybrids");
    if (outlier_num_max && delete_outliers. empty ())
      throw runtime_error ("-outlier_num_max requires -delete_outliers");
    if (! leaf_errors. empty () && noqual)
      throw runtime_error ("-noqual prohibits -leaf_errors");
    IMPLY (refresh_dissim, ! dissim_request. empty ());


    DistTree::printParam (cout);
    cout << "Root: " << (root_topological ? "topological" : "by length") << endl;
    cout << endl;


    unique_ptr<DistTree> tree;
    {
      const Chronometer_OnePass cop ("Initial topology");  
      tree. reset (isDirName (dataFName)
                     ? new DistTree (dataFName, true, true, /*optimize*/ new_only)
                     : input_tree. empty ()
                       ? new DistTree (            dataFName, dissimAttrName, multAttrName)
                       : new DistTree (input_tree, dataFName, dissimAttrName, multAttrName)
                  );
    }
    ASSERT (tree. get ());
    
    if (fix_discernibles)
    {
      cerr << "Fixing discernibles ..." << endl;
      tree->setDiscernibles ();
    }
      
    if (qc_on)
    {
      cerr << "QC ..." << endl;
      tree->qc (); 
    }

    tree->printInput (cout);
    cout << endl;
    
    
    
    bool predictionChanged = false;


    if (! deleteFName. empty ())
    {
      cerr << "Deleting ..." << endl;
      size_t deleted = 0;
      {
        LineInput f (deleteFName, 10000, 1);
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
          if (tree->optimizable ())
            prog (tree->absCriterion2str ());
          else
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
      cerr << "Keeping ..." << endl;
      VectorPtr<Leaf> toDelete;  toDelete. reserve (tree->name2leaf. size ());
      for (const auto& it : tree->name2leaf)
        toDelete << it. second;
      toDelete. sort ();
      ASSERT (toDelete. isUniq ());
      VectorPtr<Leaf> toKeep;  toKeep. reserve (tree->name2leaf. size ());
      {
        LineInput f (keepFName, 10000, 1);
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
          if (verbose ())
            tree->saveFile (output_tree);  

          if (! skip_len)            
          {
            const Chronometer_OnePass cop ("Initial arc lengths");

            const size_t lenArc_deleted = tree->optimizeLenArc ();
            cout << "# Nodes deleted = " << lenArc_deleted << endl;

            const size_t lenNode_deleted = tree->optimizeLenNode ();
            cout << "# Nodes deleted = " << lenNode_deleted << endl;
            
            tree->qc ();
            if (verbose ())
              tree->print (cout);  
            tree->reportErrors (cout);
            
            if (lenArc_deleted || lenNode_deleted)
              predictionChanged = true;
          }
          

          if (reinsert)
          {
            cerr << "Optimizing topology: re-insert ..." << endl;
            const Chronometer_OnePass cop ("Topology optimization: re-insert");
            tree->optimizeReinsert ();  
            predictionChanged = true;
          }
          
          if (predictionChanged)
         		tree->setDissimMult ();
          
          if (! skip_topology)
          {
          //cerr << "Optimizing topology: subgraphs ..." << endl;
            cout << "Optimizing topology: subgraphs ..." << endl;
            const Chronometer_OnePass cop ("Topology optimization: local");
            size_t iter_max = numeric_limits<size_t>::max ();
            if (subgraph_fast)
              minimize (iter_max, max<size_t> (1, (size_t) log2 ((Real) leaves) / areaRadius_std));  // PAR
            if (subgraph_iter_max)
            	minimize (iter_max, subgraph_iter_max);
            ASSERT (iter_max);
          //bool hybridDeleted = true;
            size_t iter = 0;
            for (; iter < iter_max /*|| (delete_all_hybrids && hybridDeleted)*/; iter++)
          	{
              cerr << "Iteration " << iter + 1;
          		if (iter_max < numeric_limits<size_t>::max ())
          	    cerr << " / " << iter_max; 
          	  cerr << " ..." << endl;
          		const Real absCriterion_old = tree->absCriterion;
              tree->optimizeDissimCoeffs ();  
              if (reinsert_iter)
                tree->optimizeReinsert ();  
              tree->optimizeLargeSubgraphs ();  
            //hybridDeleted = false;
              if (hybridF. get ())
              	/*hybridDeleted =*/ deleteHybrids (*tree, false/*??*/, hybridParentPairsF. get (), *hybridF, dissim_request. empty () ? nullptr : & hybridDissimRequests);
              tree->saveFile (output_tree); 
            #if 0
              if (hybridDeleted)
              	continue;
            #endif
              if (! absCriterion_old)
              	break;
              const Real relImprovement = (absCriterion_old - tree->absCriterion) / absCriterion_old;
              if (relImprovement <= 1e-4 / (Real) tree->name2leaf. size ())  // PAR
              	break;
              tree->setDissimMult ();
              cerr << tree->absCriterion2str () /*<< ' ' << relImprovement*/ << endl; 
            }
            cout << "# Iterations of subgraph optimization: " << iter + 1 << endl;
            tree->reportErrors (cout);
          }
          
          cerr << "Re-rooting ..." << endl;
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


      // Outliers 
      if (! delete_outliers. empty ())
      {
        cerr << "Finding criterion outliers ..." << endl;
	      tree->setLeafAbsCriterion (); 
	      Real outlier_min = NaN;
	      const VectorPtr<Leaf> outliers (tree->findCriterionOutliers (1e-10, outlier_min));  // PAR  
	      const ONumber on (cout, absCriterionDecimals, false);  
	      cout << "# Outliers: " << outliers. size () << endl;
        OFStream f (delete_outliers);
	      if (! outliers. empty ())
	      {
  	      if (verbose ())
  	      {
    	      cout << "Min. " << outlierCriterion << " of outliers: " << outlier_min << endl;
    	      for (const Leaf* leaf : outliers)
    	        cout         << leaf->name 
    	             << '\t' << leaf->getRelCriterion ()
    	             << '\t' << leaf->absCriterion  
    	             << endl;
    	    }
          cerr << "Deleting outliers ..." << endl;
          {
            Progress prog (outliers. size ());
            size_t removed = 0;
            for (const Leaf* leaf : outliers)
            {
              if (outlier_num_max && removed >= outlier_num_max)
                break;
              f << leaf->name << endl;
              tree->removeLeaf (var_cast (leaf), true);
              prog (tree->absCriterion2str ());
              removed++;
            }
            if (removed)
              tree->setDissimMult ();
          }
          tree->reportErrors (cout);
        }
        cout << endl;
      }
      tree->qc ();
    
      
      if (! noqual)
      {
        tree->setNodeAbsCriterion ();  // --> setLeafAbsCriterion() ??
        tree->qc ();
      }


      cout << "OUTPUT:" << endl;  
      tree->reportErrors (cout);
      if (verbose ())
        tree->printAbsCriterion_halves ();  
      cout << endl;
    }
    
    
    chron_getBestChange.    print (cout);
    chron_tree2subgraph.    print (cout);
    chron_subgraphOptimize. print (cout); 
    chron_subgraph2tree.    print (cout); 


    // Tree model is fixed
    

    if (reroot)
    {
      cerr << "Re-rooting ..." << endl;
      cout << "Ave. radius: " << tree->reroot (root_topological) << endl;
      cout << endl;
      tree->qc ();
    }
    else 
      if (! reroot_at. empty ())
      {
			  Tree::LcaBuffer buf;
        const DTNode* underRoot = tree->lcaName2node (reroot_at, buf);
        tree->reroot (var_cast (underRoot), underRoot->len / 2);
        if (tree->optimizable () && ! noqual)
	        tree->setNodeAbsCriterion ();  // --> setLeafAbsCriterion() ??
        tree->qc ();
      }
      

    if (! noqual && tree->optimizable ())
    {
		  const ONumber on (cout, 2, false);  // PAR
    	const Real dissim_var = tree->setErrorDensities ();  
      cout << "Relative epsilon2_0 = " << sqrt (dissim_var / tree->target2_sum) * 100 << " %" << endl;
        // Must be << "Average arc error"
      cout << "Mean residual = " << tree->getMeanResidual () << endl;
      cout << "Correlation between residual^2 and dissimilarity = " << tree->getSqrResidualCorr () << endl;  // ??
      cout << endl;
      tree->qc ();
    }
    
    
    tree->saveFile (output_tree); 
    tree->saveDissimCoeffs (output_dissim_coeff);
    tree->saveFeatureTree (output_feature_tree);


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
      checkOptimizable (*tree, "leaf_errors");
    //tree->setNodeAbsCriterion ();   // Done above
      const Dataset ds (tree->getLeafErrorDataset ());
      OFStream f (leaf_errors + dmSuff);
	    ds. saveText (f);    
    }


    if (! output_dissim. empty ())
    {
      checkOptimizable (*tree, "output_dissim");
      OFStream f (output_dissim);
      tree->saveDissim (f, output_dist_etc);
    }
    
    
    if (! dissim_request. empty ())
    {
      cerr << "Finding missing leaf pairs ..." << endl;      
      Vector<Pair<const Leaf*>> pairs (tree->getMissingLeafPairs_ancestors (sparsingDepth, refresh_dissim));
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


    if (! delete_outliers. empty ())  // Parameter is performed above
      checkOptimizable (*tree, "delete_outliers");  
    if (! delete_hybrids. empty ())  // Parameter is performed above
      checkOptimizable (*tree, "delete_hybrids");  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


