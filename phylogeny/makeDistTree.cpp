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



const string outlierCriterion ("rel. leaf error");



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Optimize topology and compute statistics for a distance tree")
	{
		// Input
	  addKey ("input_tree", "Directory with a tree of " + dmSuff + "-files ending with '/' or a tree file. If empty then neighbor-joining");
	  addKey ("data", dmSuff + "-file without \"" + dmSuff + "\", may contain more or less objects than <input_tree> does; or directory with data");
	  addKey ("dissim", "Dissimilarity attribute name in the <data> file");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
	  
	  // Processing
	  addKey ("remove", "Remove leaves whose names are in the indicated file");
	  addFlag ("sparse", "Make the initial dissimilarity matrix sparse");
	  addFlag ("topology", "Optimize topology, arc lengths and re-root");
	  addFlag ("whole", "Optimize whole topology, otherwise by subgraphs of radius " + toString (areaRadius_std));
	  addFlag ("reroot", "Re-root");
	  addFlag ("root_topological", "Root minimizes average topologcal depth, otherwise average length to leaves weighted by subtree length");
	  addKey  ("reroot_at", string ("Interior node denoted as \'A") + DistTree::objNameSeparator + "B\', which is the LCA of A and B. Re-root above the LCA in the middle of the arc");
	//addFlag ("strong_outliers", "Use log(" + outlierCriterion + ") as the outlier criterion to produce fewer of them"); 
	  addKey ("remove_outliers", "Remove outliers by " + outlierCriterion + " and save them in the indicated file");

    // Output
	  addKey ("output_tree", "Resulting tree");
	  addKey ("output_feature_tree", "Resulting tree in feature tree format");
	  addKey ("leaf_errors", "File with relative errors of leaves");
	//addKey ("pair_residuals", dmSuff + "-file with quality statistics for each object pair"); ??
	  addKey ("arc_length_stat", "File with arc length statistics: " + Tree::printArcLengthsColumns ());
	  addKey ("output_dissim", "File with dissimilarities used in the tree, tab-delimited line format: <obj1> <obj2> <dissim>");
	  addKey ("dissim_request", "File with requests to comoute dissimilarities, tab-delimited line format: <obj1> <obj2>");
	}
	
	
	
	void body () const
  {
	  const string input_tree          = getArg ("input_tree");
	  const string dataFName           = getArg ("data");
	  const string dissimAttrName      = getArg ("dissim");
	               varianceType        = str2varianceType (getArg ("variance"));  // Global    
		const string removeFName         = getArg ("remove");
		      bool   sparse              = getFlag ("sparse");
		const bool   topology            = getFlag ("topology");
		const bool   whole               = getFlag ("whole");
		const bool   reroot              = getFlag ("reroot");
		const bool   root_topological    = getFlag ("root_topological");
		const string reroot_at           = getArg ("reroot_at");
  //const bool   strong_outliers     = getFlag ("strong_outliers");
		const string remove_outliers     = getArg ("remove_outliers");
		const string output_tree         = getArg ("output_tree");
		const string output_feature_tree = getArg ("output_feature_tree");
		const string leaf_errors         = getArg ("leaf_errors");
	//const string pair_residuals      = getArg ("pair_residuals");
		const string arc_length_stat     = getArg ("arc_length_stat");
		const string output_dissim       = getArg ("output_dissim");
		const string dissim_request      = getArg ("dissim_request");
		
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
    IMPLY (whole, topology);


    DistTree::printParam (cout);
    if (sparse)
      cout << "Sparsing depth = " << sparsingDepth << endl;
    if (topology)
      cout << "Topology optimization: " << (whole ? "whole" : "subgraphs") << endl;
    cout << "Root: " << (root_topological ? "topological" : "by length") << endl;
    cout << endl;


    Common_sp::AutoPtr<DistTree> tree;
    {
      const Chronometer_OnePass cop ("Initial topology");  
      tree = isDirName (dataFName)
               ? new DistTree (dataFName, true, true)
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
    
    chron_getBestChange.    print (cout);
    chron_tree2subgraph.    print (cout);
    chron_subgraphOptimize. print (cout); 
    chron_subgraph2tree.    print (cout); 

    
    if (tree->optimizable ())
    {
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
            prog (real2str (tree->absCriterion, criterionDecimals));  
          }
        }
        tree->reportErrors (cout);
        tree->qc ();
      }


      if (topology)
      {
        const size_t leaves = tree->root->getLeavesSize ();
        if (leaves > 3)
        {
          if (verbose ())
            tree->saveFile (output_tree);  
            
          {
            const Chronometer_OnePass cop ("Initial arc lengths");

            tree->optimizeLenArc ();
            cout << "# Nodes deleted = " << tree->finishChanges () << endl;
            
            tree->optimizeLenNode ();  
            cout << "# Nodes deleted = " << tree->finishChanges () << endl;
            
            tree->qc ();
            if (verbose ())
              tree->print (cout);  
            tree->reportErrors (cout);
          }
          
          {
            cout << "Optimizing topology ..." << endl;
            const Chronometer_OnePass cop ("Topology and arc length optimization");
            if (whole)
              tree->optimizeIter (0, output_tree);
            else
            {
              tree->optimizeSubgraphs (areaRadius_std);  
                // optimizeSubtreesIter () almost does not improve
            }
            tree->reportErrors (cout);
          }
          
          cout << endl << "Ave. radius: " << tree->reroot (root_topological) << endl;
        }
        else if (leaves == 3)
          tree->optimize3 ();
        else if (leaves == 2)
          tree->optimize2 ();
      }


      // Outliers
      Real outlier_min = NAN;
      tree->setLeafAbsCriterion ();
      const VectorPtr<Leaf> outliers (tree->findOutliers (true /*strong_outliers*/, outlier_min));
      cout << endl << "# Outliers: " << outliers. size () << endl;
      cout << "Min." << (true /*strong_outliers*/ ? " log" : "") << " " << outlierCriterion << " of outliers: " << outlier_min << endl;
      for (const Leaf* leaf : outliers)
        cout         << leaf->name 
             << '\t' << leaf->getRelCriterion (true /*strong_outliers*/)
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
            prog (real2str (tree->absCriterion, criterionDecimals));  
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

      cout << "Relative epsilon2_0 = " << sqrt (tree->setErrorDensities () / tree->dissim2_sum) * 100 << " %" << endl;
        // Must be << "Average arc error"
      cout << "Mean residual = " << tree->getMeanResidual () << endl;
      cout << "Correlation between residual^2 and dissimilarity = " << tree->getSqrResidualCorr () << endl;  // ??
    }
    

    if (reroot)
      cout << endl << "Ave. radius: " << tree->reroot (root_topological) << endl;
    else 
      if (! reroot_at. empty ())
      {
        const DTNode* underRoot = tree->lcaName2node (reroot_at);
        tree->reroot (const_cast <DTNode*> (underRoot), underRoot->len / 2);
      }
      
  //tree->sort ();
    tree->setFrequentChild (rareProb);  
    tree->setFrequentDegree (rareProb); 

    tree->saveFile (output_tree);  // make size to be O(n) - extended Newick !??
    tree->saveFeatureTree (output_feature_tree);

    
    {
      const ONumber on (cout, criterionDecimals, false);
      cout << endl;
      cout << "# Interior nodes (with root) = " << tree->countInteriorNodes () << " (max = " << tree->getDiscernables (). size () - 1 << ')' << endl;
      cout << "# Interior undirected arcs = " << tree->countInteriorUndirectedArcs () << endl;
      cout << "Tree length = " << tree->getLength () << endl;
      cout << "Min. discernable leaf length = " << tree->getMinLeafLen () << endl;
        // = 0 => epsilon2_0 > 0
    #if 0
      if (sparse) 
      {
        const size_t missing = tree->selectPairs (). size ();
        cout << "Missing dissimilarities = " << missing << " (" << (Real) missing / (Real) tree->dissimSize_max () * 100 << " %)" << endl;
      }
    #endif
      cout << "Ave. arc length = " << tree->getAveArcLength () << endl;
        // Check exponential distribution ??
      cout << "Interior height = " << tree->getInteriorHeight () << endl;
      const Real bifurcatingInteriorBranching = tree->getBifurcatingInteriorBranching ();
      cout << "Bifurcating interior branching = " << bifurcatingInteriorBranching << endl;
      if (sparse) 
        cout << "# Sparsing leaves = " << pow (bifurcatingInteriorBranching, sparsingDepth + 1) << endl;
      // #dissimilarities = 2 #discernables log_2(#discernables) #sparsing_leaves

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
        f << it. first << '\t' << it. second->getRelCriterion (true /*strong_outliers*/) << endl;  
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
      OFStream f (dissim_request);
      const Set<string> pairs (tree->selectPairs ());
      for (const string& s : pairs)
      {
        string s1 (s);
        replace (s1, Tree::objNameSeparator, '\t');
        f << s1 << endl;
      }
    }

    if (! remove_outliers. empty ())  // Parameter is performed above
      checkOptimizable (*tree, "remove_outliers");  
    if (! removeFName. empty ())  // Parameter is performed above  
      checkOptimizable (*tree, "remove");  // ??
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


