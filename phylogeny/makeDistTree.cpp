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



const string outlierCriterion ("rel. average leaf error");



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Optimize topology and compute statistics for a least-squares distance tree")
	{
		// Input
	  addKey ("input_tree", "Directory with a tree of " + dmSuff + "-files ending with '/' or a tree file. If empty then neighbor-joining");
	  addKey ("data", dmSuff + "-file without \"" + dmSuff + "\", may contain more or less objects than <input_tree> does; or directory with data for an incremental tree");
	  addKey ("dissim", "Dissimilarity attribute name in the <data> file");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
	  addKey ("dist_request", "File with requests to compute tree distances, tab-delimited line format: <obj1> <obj2>");
	  
	  // Processing
	  addKey ("remove", "Remove leaves whose names are in the indicated file");
	  addFlag ("sparse", "Make the initial dissimilarity matrix sparse");

	  addFlag ("optimize", "Optimize topology, arc lengths and re-root");
	  addFlag ("whole", "Optimize whole topology, otherwise by subgraphs of radius " + toString (areaRadius_std));
	  addFlag ("subgraph_fast", "Limit the number of iterations over the whole tree for subgraph optimizations");
	  addFlag ("skip_len", "Skip length-only optimization");
	  addFlag ("reinsert", "Re-insert subtrees");
	  addFlag ("skip_topology", "Skip topology optimization");
	  
	  addFlag ("new_only", "Optimize only new objects in an incremental tree, implies not -optimize");

	  addFlag ("reroot", "Re-root");
	  addFlag ("root_topological", "Root minimizes average topologcal depth, otherwise average length to leaves weighted by subtree length");
	  addKey  ("reroot_at", string ("Interior node denoted as \'A") + DistTree::objNameSeparator + "B\', which is the LCA of A and B. Re-root above the LCA in the middle of the arc");
	  addKey ("remove_outliers", "Remove outliers by " + outlierCriterion + " and save them in the indicated file");

    // Output
	  addKey ("output_tree", "Resulting tree");
	  addKey ("output_feature_tree", "Resulting tree in feature tree format");
	  addKey ("leaf_errors", "File with relative errors of leaves");
	//addKey ("pair_residuals", dmSuff + "-file with quality statistics for each object pair"); ??
	  addKey ("arc_length_stat", "File with arc length statistics: " + Tree::printArcLengthsColumns ());
	  addKey ("output_dissim", "File with dissimilarities used in the tree, tab-delimited line format: <obj1> <obj2> <dissim>");
	  addKey ("dissim_request", "File with requests to compute needed dissimilarities, tab-delimited line format: <obj1> <obj2>");
	  addKey ("output_dist", "File with tree distances, tab-delimited line format: <obj1> <obj2> <dist>");
	}
	
	
	
	void body () const final
  {
	  const string input_tree          = getArg ("input_tree");
	  const string dataFName           = getArg ("data");
	  const string dissimAttrName      = getArg ("dissim");
	               varianceType        = str2varianceType (getArg ("variance"));  // Global
		const string dist_request        = getArg ("dist_request");
	               
		const string removeFName         = getArg ("remove");
		      bool   sparse              = getFlag ("sparse");
		      
		const bool   optimize            = getFlag ("optimize");
		const bool   whole               = getFlag ("whole");
		const bool   subgraph_fast       = getFlag ("subgraph_fast");
		const bool   skip_len            = getFlag ("skip_len");
		const bool   reinsert            = getFlag ("reinsert");
		const bool   skip_topology       = getFlag ("skip_topology");
		const bool   new_only            = getFlag ("new_only");
		
		const bool   reroot              = getFlag ("reroot");		
		const bool   root_topological    = getFlag ("root_topological");
		const string reroot_at           = getArg ("reroot_at");
		const string remove_outliers     = getArg ("remove_outliers");
		
		const string output_tree         = getArg ("output_tree");
		const string output_feature_tree = getArg ("output_feature_tree");
		const string leaf_errors         = getArg ("leaf_errors");
	//const string pair_residuals      = getArg ("pair_residuals");
		const string arc_length_stat     = getArg ("arc_length_stat");
		const string output_dissim       = getArg ("output_dissim");
		const string dissim_request      = getArg ("dissim_request");
		const string output_dist         = getArg ("output_dist");
		
	  IMPLY (isDirName (input_tree), ! sparse);
		ASSERT (! (reroot && ! reroot_at. empty ()));
    if (isDirName (dataFName))
    {
      if (! input_tree. empty ())
        throw runtime_error ("Input tree must be in " + dataFName);
      if (! dissimAttrName. empty ())
        throw runtime_error ("Non-empty dissimilarity attribute with no " + dmSuff + "-file");
      if (sparse)
        throw runtime_error ("Further sparsing of " + dataFName + " cannot be done");
      sparse = true;
    }
    else
      if (dataFName. empty () != dissimAttrName. empty ())
        throw runtime_error ("The both data file and the dissimilarity attribute must be either present or absent");
    if (! dist_request. empty () && output_dist. empty ())
    	throw runtime_error ("dist_request exist, but no output_dist");
    IMPLY (whole,         optimize);
    IMPLY (subgraph_fast, optimize);
    IMPLY (skip_len,      optimize);
    IMPLY (reinsert,      optimize);
    IMPLY (skip_topology, optimize);
    IMPLY (subgraph_fast, ! whole);
    IMPLY (new_only, ! optimize);


    DistTree::printParam (cout);
    if (sparse)
      cout << "Sparsing depth = " << sparsingDepth << endl;
    cout << "Root: " << (root_topological ? "topological" : "by length") << endl;
    cout << endl;


    Common_sp::AutoPtr<DistTree> tree;
    {
      const Chronometer_OnePass cop ("Initial topology");  
      tree = isDirName (dataFName)
               ? new DistTree (dataFName, true, true, new_only)
               : input_tree. empty ()
                 ? new DistTree (dataFName, dissimAttrName, sparse)
                 : isDirName (input_tree)
                   ? new DistTree (input_tree, dataFName, dissimAttrName)
                   : new DistTree (input_tree, dataFName, dissimAttrName, sparse);
    }
    ASSERT (tree. get ());
    tree->qc ();     

    tree->printInput (cout);
    cout << endl;


    if (! removeFName. empty ())
    {
      cout << endl << "Removing ..." << endl;
      {
        LineInput f (removeFName, 10000, 1);
        Progress prog;
        while (f. nextLine ())
        {
          trim (f. line);
          const Leaf * leaf = findPtr (tree->name2leaf, f. line);
          if (! leaf)
          {
            cout << "Leaf " << f. line << " not found" << endl;
            continue;
          }
          tree->removeLeaf (const_cast <Leaf*> (leaf));
          if (tree->optimizable ())
            prog (tree->absCriterion2str ());
          else
          	prog ();
        }
      }
	    if (tree->optimizable ())
        tree->reportErrors (cout);
      tree->qc ();
    }


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
          }
          
          if (reinsert)
          {
            cout << "Optimizing topology: re-insert ..." << endl;
            const Chronometer_OnePass cop ("Topology optimization: re-insert");
            tree->optimizeReinsert ();  
          }
          
          if (! skip_topology)
          {
            cout << string ("Optimizing topology: ") + (whole ? "neighbors" : "subgraphs") + " ..." << endl;
            const Chronometer_OnePass cop ("Topology optimization: local");
            if (whole)
              tree->optimizeWholeIter (0, output_tree);
            else
            {
              const size_t iter_max = (size_t) log2 ((Real) leaves) / areaRadius_std + 1;  // PAR
              size_t iter = 0;
            	for (;;)
            	{
	              iter++;
            		if (subgraph_fast)
            		{
            			if (iter > iter_max)
            			  break;
            	    cout << "Iteration " << iter << '/' << iter_max << " ..." << endl;
            	  }
            		const Real absCriterion_old = tree->absCriterion;
	              tree->optimizeSubgraphs (areaRadius_std);  
	              if (tree->absCriterion == 0)
	              	break;
	              if ((absCriterion_old - tree->absCriterion) / tree->absCriterion < 1e-6)  // PAR
	              	break;
	            }
            }
            tree->reportErrors (cout);
          }
          
          const Real radius_ave = tree->reroot (root_topological);
          cout << endl << "Ave. radius: " << radius_ave << endl;
        }
        else if (leaves == 3)
          tree->optimize3 ();
        else if (leaves == 2)
          tree->optimize2 ();
      }


      // Outliers
      Real outlier_min = NAN;
      tree->setLeafAbsCriterion ();
      const VectorPtr<Leaf> outliers (tree->findCriterionOutliers (0.1, outlier_min));  // PAR
      cout << endl << "# Outliers: " << outliers. size () << endl;
      cout << "Min. " << outlierCriterion << " of outliers: " << outlier_min << endl;
      for (const Leaf* leaf : outliers)
        cout         << leaf->name 
             << '\t' << leaf->getRelCriterion ()
             << '\t' << leaf->absCriterion
             << endl;
      if (! remove_outliers. empty ())
      {
        cout << endl << "Removing outliers ..." << endl;
        {
          OFStream f (remove_outliers);
          Progress prog ((uint) outliers. size ());
          for (const Leaf* leaf : outliers)
          {
            f << leaf->name << endl;
            tree->removeLeaf (const_cast <Leaf*> (leaf));
            prog (tree->absCriterion2str ());
          }
        }
        tree->reportErrors (cout);
        tree->qc ();
      }

      
      cout << endl << "OUTPUT:" << endl;  
      tree->reportErrors (cout);
      tree->printAbsCriterion_halves ();  // skip if isDirName(dataFName) ??
      tree->setLeafAbsCriterion ();
      tree->qc ();
    }
    

    if (reroot)
      cout << endl << "Ave. radius: " << tree->reroot (root_topological) << endl;
    else 
      if (! reroot_at. empty ())
      {
        const DTNode* underRoot = tree->lcaName2node (reroot_at);
        tree->reroot (const_cast <DTNode*> (underRoot), underRoot->len / 2);
      }
      

    if (tree->optimizable ())
    {
    	const Real dissim_var = tree->setErrorDensities ();
      cout << "Relative epsilon2_0 = " << sqrt (dissim_var / tree->dissim2_sum) * 100 << " %" << endl;
        // Must be << "Average arc error"
      cout << "Mean residual = " << tree->getMeanResidual () << endl;
      cout << "Correlation between residual^2 and dissimilarity = " << tree->getSqrResidualCorr () << endl;  // ??
    }
    tree->qc ();
    
    
    chron_getBestChange.    print (cout);
    chron_tree2subgraph.    print (cout);
    chron_subgraphOptimize. print (cout); 
    chron_subgraph2tree.    print (cout); 


  //tree->sort ();
    tree->setFrequentChild (rareProb);  
    tree->setFrequentDegree (rareProb); 

    tree->saveFile (output_tree);  
    tree->saveFeatureTree (output_feature_tree);


  #if 0
    {
      Real arcLen_min = NAN;
      Real outlier_EValue_max = 10;  // ??
      while (outlier_EValue_max >= 1e-6)
      {
        const VectorPtr<DTNode> tooLongArcs (tree->findOutlierArcs (outlier_EValue_max, arcLen_min));
        cout << endl;
        {
          ONumber on (cout, 1, true);
          cout << "outlier_EValue_max: " << outlier_EValue_max << endl;
        }
        cout << "# Too long arcs: " << tooLongArcs. size () << endl;
        cout << "Min. length of too long arcs: " << arcLen_min << endl;
        outlier_EValue_max /= 10;
      }
    }
  #endif
  //tree->findTopologicalClusters ();
  #if 0
    {
      Real outlier_EValue_max = 0.001;  
      while (outlier_EValue_max >= 1e-10)
      {
        ONumber on (cout, 1, true);
        cout << "outlier_EValue_max: " << outlier_EValue_max << endl;
        tree->findDepthOutliers (outlier_EValue_max);  
        outlier_EValue_max /= 10;
      }
    }
  #endif

    
    // Statistics
    {
      const ONumber on (cout, criterionDecimals, false);
      cout << endl;
      cout << "# Interior nodes (with root) = " << tree->countInteriorNodes () << " (max = " << tree->getDiscernibles (). size () - 1 << ')' << endl;
      cout << "# Interior undirected arcs = " << tree->countInteriorUndirectedArcs () << endl;
      cout << "Tree length = " << tree->getLength () << endl;
      {
      	const ONumber on1 (cout, 3, true);  // PAR
        cout << "Min. discernible leaf length = " << tree->getMinLeafLen () << endl;
          // = 0 => epsilon2_0 > 0
      }
      cout << "Ave. arc length = " << tree->getAveArcLength () << endl;
        // Check exponential distribution ??
      cout << "Interior height = " << tree->getInteriorHeight () << endl;
      const Real bifurcatingInteriorBranching = tree->getBifurcatingInteriorBranching ();
      cout << "Bifurcating interior branching = " << bifurcatingInteriorBranching << endl;
      if (sparse) 
        cout << "# Sparsing leaves = " << pow (bifurcatingInteriorBranching, sparsingDepth + 1) << endl;
      // #dissimilarities = 2 #discernibles log_2(#discernibles) #sparsing_leaves

      {      
        size_t freqChildrenInteriors = 0;
        size_t freqChildrenLeaves    = 0;
        size_t stableInteriors       = 0;
        size_t stableLeaves          = 0;
        for (const DiGraph::Node* node : tree->nodes)
        {
          const Tree::TreeNode* tn = static_cast <const Tree::TreeNode*> (node);
          if (tn->frequentChild)
          {
            if (tn->isInteriorType ())
              freqChildrenInteriors++;
            else
              if (tn->isLeafType ())
                freqChildrenLeaves++;
          }
          if (   tn->isInteriorType ()
              && tn->frequentDegree >= 3
             )
            stableInteriors++;
          if (   tn->isLeafType ()
              && tn->frequentDegree == 1
             )
            stableLeaves++;
        }
        
        // Depend on root ??
        cout << "# Frequent children interior nodes = " << freqChildrenInteriors << endl;
        cout << "# Frequent children leaves = "         << freqChildrenLeaves << endl;
        
        cout << "# Frequent interior nodes = "          << stableInteriors << endl;
        cout << "# Frequent leaves = "                  << stableLeaves << endl;
        cout << "Rareness threshold = " << rareProb * 100 << " %" << endl;
      }      
    }

    
    if (! leaf_errors. empty ())
    {
      checkOptimizable (*tree, "leaf_errors");
      OFStream f (leaf_errors);
    //tree->setLeafAbsCriterion ();   // Done above
      for (const auto& it : tree->name2leaf)
        f << it. first << '\t' << it. second->getRelCriterion () << endl;  
    }

  #if 0
    ??
    if (! pair_residuals. empty ())
    {
      checkOptimizable (*tree, "pair_residuals");
      OFStream f ("", pair_residuals, dmExt);
      const RealAttr1* resid2Attr  = tree->getResiduals2 ();
      const RealAttr1* logDiffAttr = tree->getLogPredictionDiff ();
      tree->pairResiduals2dm (resid2Attr, logDiffAttr, f); 
    }
  #endif
    
    if (! arc_length_stat. empty ())
    {
      // dm-file ??
      // cout << arcLenRel.SD after outlier removing ??
      OFStream f (arc_length_stat);
      const ONumber on (f, dissimDecimals, true);
      tree->printArcLengths (f);
    }
    
    if (! output_dissim. empty ())
    {
      checkOptimizable (*tree, "output_dissim");
      OFStream f (output_dissim);
      tree->saveDissim (f);
    }
    
    
    if (! dissim_request. empty ())
    {
      cout << endl << "Finding missing leaf pairs ..." << endl;
      tree->setReprLeaves ();
      tree->dissims. sort ();
      
      Vector<Pair<const Leaf*>> pairs (tree->getMissingLeafPairs_ancestors (sparsingDepth));
      cout << "# Ancestor-based requests: " << pairs. size () << endl;

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
      
      pairs. sort ();
      pairs. uniq ();      
      
      cout << "# Total requests: " << pairs. size () << endl;
    #endif
      {
        OFStream f (dissim_request);      
        for (const auto& p : pairs)
          f << p. first->name << '\t' << p. second->name << endl;
      }
    }


    Vector<Pair<const Leaf*>> distRequestPairs;  distRequestPairs. reserve (tree->dissims. size ());
    if (! dist_request. empty ())
    {
    	LineInput f (dist_request);
    	while (f. nextLine ())
    	{
    		istringstream iss (f. line);
    		string name1, name2;
    		iss >> name1 >> name2;
        const Leaf* leaf1 = findPtr (tree->name2leaf, name1);
        if (! leaf1)
        	throw runtime_error ("Object " + name1 + " is not found");
        const Leaf* leaf2 = findPtr (tree->name2leaf, name2);
        if (! leaf2)
        	throw runtime_error ("Object " + name2 + " is not found");
			  distRequestPairs << Pair<const Leaf*> (leaf1, leaf2);
    	}
    }
    
    
    if (! output_dist. empty ())
    {
    	// distRequestPairs
    	if (dist_request. empty ())
    	{
	    	VectorPtr<Leaf> leaves;  leaves. reserve (tree->nodes. size ());
	 	    for (const DiGraph::Node* node : tree->nodes)
	        if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
	        	leaves << leaf;
	      FFOR (size_t, i, leaves. size ())
	        FOR (size_t, j, i)
	        {
	        	const Leaf* leaf1 = leaves [i];
	        	const Leaf* leaf2 = leaves [j];
					  if (leaf1->name > leaf2->name)
					    swap (leaf1, leaf2);
					  distRequestPairs << Pair<const Leaf*> (leaf1, leaf2);
					}
    	}

    	OFStream f (output_dist);
    	for (const auto p : distRequestPairs)
    	{
    		const Leaf* leaf1 = p. first;
    		const Leaf* leaf2 = p. second;
			  const Tree::TreeNode* lca_ = nullptr;
			  const VectorPtr<Tree::TreeNode> path (Tree::getPath (leaf1, leaf2, nullptr, lca_));
		  	f         << leaf1->name 
		  	  << '\t' << leaf2->name
		  	  << '\t' << DistTree::path2prediction (path)
		  	  << endl;
		  }
    }


    if (! remove_outliers. empty ())  // Parameter is performed above
      checkOptimizable (*tree, "remove_outliers");  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


