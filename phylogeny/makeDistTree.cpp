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
	  addFlag ("sparse", "Make the initial dissimilarity matrix sparse");
	  addFlag ("topology", "Optimize topology, arc lengths and re-root");
	  addFlag ("whole", "Optimize whole topology, otherwise by subtrees of radius " + toString (areaRadius_std));
	  addFlag ("reroot", "Re-root");
	  addKey  ("reroot_at", string ("Interior node denoted as \'A") + DistTree::objNameSeparator + "B\', which is the LCA of A and B. Re-root above the LCA in the middle of the arc");
	//addFlag ("sparse_add", "Add the dissimilarities to the dissimilarity matrix sparsely");  

    // Output
	  addKey ("output_tree", "Resulting tree");

	  addKey ("output_feature_tree", "Resulting tree in feature tree format");
	  addKey ("leaf_errors", "File with relative errors of leaves");
	  addKey ("pair_residuals", dmSuff + "-file with quality statistics for each object pair");
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
		const bool topology              = getFlag ("topology");
		const bool whole                 = getFlag ("whole");
		const bool reroot                = getFlag ("reroot");
		const string reroot_at           = getArg ("reroot_at");
		const bool sparse                = getFlag ("sparse");
	//const bool sparse_add            = getFlag ("sparse_add");
		const string output_tree         = getArg ("output_tree");
		const string output_feature_tree = getArg ("output_feature_tree");
		const string leaf_errors         = getArg ("leaf_errors");
		const string pair_residuals      = getArg ("pair_residuals");
		const string arc_length_stat     = getArg ("arc_length_stat");
		const string output_dissim       = getArg ("output_dissim");
		const string dissim_request      = getArg ("dissim_request");
		
	  IMPLY (isRight (input_tree, "/"), ! sparse);
		ASSERT (! (reroot && ! reroot_at. empty ()));
    if (isRight (dataFName, "/"))
    {
      if (! input_tree. empty ())
        throw runtime_error ("Input tree must be in " + dataFName);
      if (! dissimAttrName. empty ())
        throw runtime_error ("Non-empty dissimilarity attribute with no " + dmSuff + "-file");
      if (sparse)
        throw runtime_error ("Sparsing " + dataFName + " cannot be done");
    }
    else
      if (dataFName. empty () != dissimAttrName. empty ())
        throw runtime_error ("The both data file and the dissimilarity attribute must be either present or absent");


    DistTree::printParam (cout);
    if (topology)
    {
      cout << "Topology optimization: " << (whole ? "whole" : "subgraphs") << endl;
      if (sparse)
        cout << "Sparsing depth = " << sparsingDepth << endl;
    }
    cout << endl;


    Common_sp::AutoPtr<DistTree> tree;
    {
      Chronometer_OnePass cop ("Initial topology");
      tree = isRight (dataFName, "/")
               ? new DistTree (dataFName, true)
               : input_tree. empty ()
                 ? new DistTree (dataFName, dissimAttrName, sparse)
                 : isRight (input_tree, "/")
                   ? new DistTree (input_tree, dataFName, dissimAttrName)
                   : new DistTree (input_tree, dataFName, dissimAttrName, sparse);
    }
    ASSERT (tree. get ());
  //if (verbose ())
      tree->qc ();     

    tree->printInput (cout);
    cout << endl;

    
    if (tree->optimizable ())
    {
      if (topology)
      {
        const size_t leaves = tree->root->getLeavesSize ();
        if (leaves > 3)
        {
          if (verbose ())
            tree->saveFile (output_tree);  
            
          if (! isRight (dataFName, "/"))
          {
            Chronometer_OnePass cop ("Initial arc lengths");

          #if 0
            tree->reportErrors (cout);  
            cout << endl;
          #endif

          #if 0
            EXEC_ASSERT (tree->optimizeLenAll ());
            tree->reportErrors (cout);  
            cout << "# Nodes deleted = " << tree->finishChanges () << endl;
            cout << endl;
          #endif

            EXEC_ASSERT (tree->optimizeLenArc ());
            cout << "# Nodes deleted = " << tree->finishChanges () << endl;
            cout << endl;
            
            tree->optimizeLenNode ();  
            cout << "# Nodes deleted = " << tree->finishChanges () << endl;
            
            tree->qc ();
            if (verbose ())
              tree->print (cout);  
            tree->reportErrors (cout);
          }
          
          {
            cout << "Optimizing topology ..." << endl;
            Chronometer_OnePass cop ("Topology and arc length optimization");
            if (whole)
              tree->optimizeIter (output_tree);
            else
              tree->optimizeSubgraphs ();  
                // optimizeSubtreesIter () almost does not improve
          }
          
          tree->reroot ();  
        }
        else if (leaves == 3)
          tree->optimize3 ();
        else if (leaves == 2)
          tree->optimize2 ();
      }

      
    #if 0
      if (tree->dissimAttr)  
      {
        tree->optimizeAdd (sparse_add, output_tree);  
        tree->reroot ();  
      //if (verbose ())
          tree->qc ();
      }
      // tree and dist-matrix match
    #endif

      if (reroot)
        tree->reroot ();
      if (! reroot_at. empty ())
      {
        const DTNode* underRoot = tree->lcaName2node (reroot_at);
        tree->reroot (const_cast <DTNode*> (underRoot), underRoot->len / 2);
      }
      
      cout << "OUTPUT:" << endl;  
      tree->reportErrors (cout);
      tree->printAbsCriterion_halves ();  
      tree->setHeight ();
      tree->setLeafAbsCriterion ();
    //if (verbose ())
        tree->qc ();

      cout << "Relative epsilon2_0 = " << sqrt (tree->setErrorDensities () / tree->dissim2_sum) * 100 << " %" << endl;
        // Must be << "Average arc error"
      cout << "Mean residual = " << tree->getMeanResidual () << endl;
      cout << "Correlation between residual^2 and dissimilarity = " << tree->getSqrResidualCorr () << endl;  // ??

      cout << endl << "Outliers:" << endl;
      const size_t outliers = tree->printLeafRelLenErros (cout, 3);  // PAR
      cout << "# Outliers: " << outliers << endl;
    }
    

    {
      Unverbose unv;
      if (verbose ())
        tree->ds. print (cout);
    }
    if (verbose ())
    {
      tree->setPrediction ();
      tree->checkAbsCriterion ("setPrediction");
    }
  
  //tree->sort ();
    tree->setFrequentChild (rareProb);  
    tree->setFrequentDegree (rareProb); 

    tree->saveFile (output_tree);
    tree->saveFeatureTree (output_feature_tree);
    
    {
      const ONumber on (cout, 4, false);
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
        cout << "# Frequent children interior nodes = " << freqChildrenInteriors << endl;
        cout << "# Frequent children leaves = "         << freqChildrenLeaves << endl;
        cout << "# Frequent interior nodes = "          << stableInteriors << endl;
        cout << "# Frequent leaves = "                  << stableLeaves << endl;
          // Incertae sedis ??
        cout << "Rareness threshold = " << rareProb * 100 << " %" << endl;
      }      
    }
      
    if (! leaf_errors. empty ())
    {
      OFStream f (leaf_errors);
      tree->setLeafAbsCriterion ();
      tree->printLeafRelLenErros (f, 0); 
    }

    if (! pair_residuals. empty ())
    {
      OFStream f ("", pair_residuals, dmExt);
      const RealAttr1* resid2Attr   = tree->getResiduals2 ();
      const RealAttr1* logDiffAttr = tree->getLogPredictionDiff ();
      tree->pairResiduals2dm (resid2Attr, logDiffAttr, f); 
    }
    
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
      OFStream f (output_dissim);
      tree->printDissim (f);
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
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


