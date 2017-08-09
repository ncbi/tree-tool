// makeDistTree.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
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
	  addPositional ("input_tree", "Directory with a tree of " + dmSuff + "-files ending with '/' or a tree file. If empty then neighbor-joining");
	  addKey ("data", dmSuff + "-file without \"" + dmSuff + "\", may contain more or less objects than <input_tree> does");
	  addKey ("dissim", "Dissimilarity attribute name in the <data> file");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
	  addFlag ("topology", "Optimize topology, arc lengths and re-root");
	  addFlag ("whole", "Optimize whole topology, otherwise by subtrees of radius " + toString (areaRadius_std));
	  addFlag ("reroot", "Re-root");
	  addKey  ("reroot_at", "Interior node denoted as \'A-B\', which is the LCA of A nd B. Re-root above the LCA in the middle of the arc");
	  addFlag ("sparse_init", "Make the initial dissimilarity matrix sparse");
	  addFlag ("sparse_add", "Add the dissimilarities to the dissimilarity matrix sparsely");

    // Output
	  addKey ("output_tree", "Resulting tree");

	  addKey ("output_feature_tree", "Resulting tree in feature tree format");
	  addKey ("leaf_errors", "File with relative errors of leaves");
	  addKey ("pair_residuals", dmSuff + "-file with quality statistics for each object pair");
	  addKey ("arc_length_stat", "File with arc length statistics: " + Tree::printArcLengthsColumns ());
	//addKey ("patr_dist", "File with patristic distances in format: <leaf name1> <leaf name2> <distance>, where <leaf name1> < <leaf name2>"); ??
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
		      string reroot_at           = getArg ("reroot_at");
		const bool sparse_init           = getFlag ("sparse_init");
		const bool sparse_add            = getFlag ("sparse_add");
		const string output_tree         = getArg ("output_tree");
		const string output_feature_tree = getArg ("output_feature_tree");
		const string leaf_errors         = getArg ("leaf_errors");
		const string pair_residuals      = getArg ("pair_residuals");
		const string arc_length_stat     = getArg ("arc_length_stat");
		IMPLY (sparse_init, ! input_tree. empty () && ! isRight (input_tree, "/"));
		ASSERT (! (reroot && ! reroot_at. empty ()));
    if (dataFName. empty () != dissimAttrName. empty ())
      throw runtime_error ("The both data file and the dissimilarity attribute must be either present or absent");


    Common_sp::AutoPtr<DistTree> tree;
    {
      Chronometer_OnePass cop ("Initial topology");
      tree = input_tree. empty ()
               ? new DistTree (dataFName, dissimAttrName)
               : isRight (input_tree, "/")
                 ? new DistTree (input_tree, dataFName, dissimAttrName)
                 : new DistTree (input_tree, dataFName, dissimAttrName, sparse_init);
    }
    ASSERT (tree. get ());
    if (verbose ())
      tree->qc ();     
      
    tree->printInput (cout);
    cout << "Max. possible dissimilarity = " << dissim_max () << endl;
    
    if (tree->optimizable ())
    {
      if (topology)
      {
        cout << "Topology optimization: " << (whole ? "whole" : "by subgraphs") << endl;
        cout << endl;
        const size_t leaves = tree->root->getLeavesSize ();
        if (leaves > 3)
        {
          if (verbose ())
            tree->saveFile (output_tree);  
          {
            Chronometer_OnePass cop ("Initial arc lengths");
            EXEC_ASSERT (tree->optimizeLen ());
            tree->finishChanges (); 
            tree->optimizeLenLocal ();  
            tree->finishChanges (); 
            if (verbose ())
            {
              tree->qc ();
              tree->print (cout);  
            }
            tree->reportErrors (cout);
          }
          {
            Chronometer_OnePass cop ("Topology and arc length optimization");
            if (whole)
              tree->optimizeIter (output_tree);
            else
              tree->optimizeSubtrees ();  
              // optimizeSubtreesIter () almost does not improve
          }
          tree->reroot ();  
        }
        else if (leaves == 3)
          tree->optimize3 ();
        else if (leaves == 2)
          tree->optimize2 ();
      }
      
      if (tree->dissimAttr)  
      {
        tree->optimizeAdd (sparse_add, output_tree);  
        tree->reroot ();  
        if (verbose ())
          tree->qc ();
      }
      // tree and dist-matrix match

      if (reroot)
        tree->reroot ();
      if (! reroot_at. empty ())
      {
        const DTNode* underRoot = tree->lcaName2node (reroot_at);
        tree->reroot (const_cast <DTNode*> (underRoot), underRoot->len / 2);
      }
        
      cout << endl;      
      tree->reportErrors (cout);
      tree->printAbsCriterion_halves ();  
      tree->setHeight ();
      tree->setLeafAbsCriterion ();
      if (verbose ())
        tree->qc ();

      cout << "Relative epsilon2_0 = " << sqrt (tree->setErrorDensities () / tree->dissim2_sum) * 100 << " %" << endl;
        // Must be << "Average arc error"
      cout << "Mean residual = " << tree->getMeanResidual () << endl;
      cout << "Correlation between residual^2 and dissimilarity = " << tree->getSqrResidualCorr () << endl;  // ??

      cout << endl << "Outliers:" << endl;
      const size_t outliers = tree->printLeafRelLenErros (cout, 3);  // PAR
      cout << "# Outliers: " << outliers << endl << endl;
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
    tree->saveFile (output_tree);
    tree->saveFeatureTree (output_feature_tree);
    
  //ONumber on (cout, 6, false);
    cout << "# Interior nodes (with root) = " << tree->size (false) << endl;
    cout << "# Interior undirected arcs = " << tree->interiorUndirectedArcs () << endl;
    cout << "Tree length = " << tree->getLength () << endl;
    cout << "Min. discernable leaf length = " << tree->getMinLeafLen () << endl;
      // = 0 => epsilon2_0 > 0
      
    if (! leaf_errors. empty ())
    {
      OFStream f ("", leaf_errors, "");
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
      OFStream f ("", arc_length_stat, "");
      ONumber on (f, 6, false);
      tree->printArcLengths (f);
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


