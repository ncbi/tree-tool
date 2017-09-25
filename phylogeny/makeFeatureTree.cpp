// makeFeatureTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "featureTree.hpp"
using namespace FeatureTree_sp;



namespace 
{



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Optimize a feature tree")
	{
		// Input
	  addKey ("input_tree", "Input file with the tree");
	  addKey ("gene_dir", "Input directory with genes for each genome");
	  addKey ("phen_dir", "Input directory with phenotypes for each genome");
	  addKey ("input_core", "Input file with root core feature ids");
	  addFlag ("no_genes", "No gene data");
	    
	  addFlag ("use_time", "Use time for MLE, otherwise parsimony method");
	//addFlag ("bad_nodes_to_root", "Move bad nodes to root");
	  addKey ("optim_iter_max", "# Iterations for tree optimization; -1: optimize time only", "0");

    // Output	    
	  addKey ("output_tree", "Output file with the tree");
	  addFlag ("set_node_ids", "Set node ids in the output tree");  
	  addKey ("output_core", "Output file with root core feature ids");
	    
	#if 0
	  addFlag ("save_phen_only", "Save phenotype changes only");
	  addKey ("save_feature", "Save changes only for a feature");
	  addKey ("save_features", "Save changes only for a set of features in the file with line format: <Feature name> <Print name>");
	#endif

	  addKey ("newick", "Output file with the tree in the Newick format");
	  addFlag ("min_newick_name", "Minimal leaf names in Newick");
	  addFlag ("qual", "Print the quality statistics measured by phenotype consistency");
	  addKey ("phen_gains", "File name to save nodes where phenotypes are gained");
	  addKey ("arc_length_stat", "File with arc length statistics in format " + Tree::printArcLengthsColumns ());
	  addKey ("patr_dist", "File with patristic distances in format: <leaf name1> <leaf name2> <distance>, where <leaf name1> < <leaf name2>");
	}



	void body () const
  {
		const string input_tree      = getArg ("input_tree");
		const string gene_dir        = getArg ("gene_dir");
		const string phen_dir        = getArg ("phen_dir");
		const string input_core      = getArg ("input_core");
		const bool no_genes          = getFlag ("no_genes");

		const bool use_time          = getFlag ("use_time");  
		const int optim_iter_max     = str2<int> (getArg ("optim_iter_max"));

		const string output_tree     = getArg ("output_tree");
		const bool set_node_ids      = getFlag ("set_node_ids");  
		const string output_core     = getArg ("output_core");

	#if 0
		const bool save_phen_only    = getFlag ("save_phen_only");  
		const string save_feature    = getArg ("save_feature");
		const string save_features   = getArg ("save_features");
  #endif

		const string newick          = getArg ("newick");
		const bool min_newick_name   = getFlag ("min_newick_name");
		const bool qual              = getFlag ("qual");  
		const string phen_gains      = getArg ("phen_gains");  
		const string arc_length_stat = getArg ("arc_length_stat");
		const string patrDistFName   = getArg ("patr_dist");

		IMPLY (! input_core. empty () && ! output_core. empty (), input_core != output_core);
	//ASSERT (save_feature. empty () || save_features. empty ());
	//ASSERT (species_id > 0);
		ASSERT (optim_iter_max >= -1);
		
		
    FeatureTree tree (/*species_id,*/ input_tree, gene_dir, phen_dir, input_core, ! no_genes);
    tree. printInput (cout);
    tree. qc ();    
    
   	if (use_time && tree. allTimeZero)
    {
   		tree. useTime (input_core);
      cout << endl;
      tree. printInput (cout);
    }
    
  #if 0
    if (bad_nodes_to_root)
    {
      const size_t n = tree. badNodesToRoot ();
      cout << endl;
      cout << "# Bad nodes moved to root: " << n << endl;
    }
  #endif
    
    if (optim_iter_max)  
    {
      if (optim_iter_max > 0)  
      {
  	  	cerr << "Optimizing ..." << endl;
        int iter = 0; 
        for (;;)
        {
    	    cout << endl;
    	    tree. dump (output_tree, true);  // Redundant
          if (iter >= optim_iter_max)
          	break;
          iter++;
    	    cout << endl << "Iter = " << iter << endl;
          if (! tree. optimize ())
          	break;	    
        }
        tree. dump (output_tree, true);  // Redundant
      }
      else
      {
  	  	cout << "Optimizing time ..." << endl;
        tree. optimizeTime ();
      }
    }

    if (! output_core. empty ())
    {
      if (tree. allTimeZero)
      {
        cout << "Finding root ..." << endl;
        tree. findRoot ();
      }
      else
      {
        cout << endl;
        cout << "Adjusting root core ..." << endl;
        tree. setCore ();
        size_t coreChange [2/*core2nonCore*/];
        tree. resetRootCore (coreChange);
        cout << "# Core to non-core: " << coreChange [true]  << endl;
        cout << "# Non-core to core: " << coreChange [false] << endl;
      }
      // Output
      tree. saveRootCore (output_core);
    }

    // Output
    tree. qc ();
  #if 0
    if (save_phen_only)
      tree. savePhenChangesOnly = true;
    if (! save_feature. empty ())
    {
      TargetFeature tf ( save_feature
	                     , "GENE"
	                     , tree. features. binSearch (Feature (save_feature, true))
	                     );
	    ASSERT (tf. index != NO_INDEX);
	    tree. targetFeatures << tf;
      tree. savePhenChangesOnly = true;
	  }
    tree. loadTargetFeatures (save_features);	  
  #endif
    tree. dump (output_tree, set_node_ids);


    if (qual)
    {    
      // Input: Feature::{genomes,gains,losses}
      cout << endl;
      cout << "Phenotype gains:" << endl;
      size_t monoPhens = 0;  // Monophyletic
      size_t phenGains = 0;
      size_t phenLosses = 0;
      size_t optionalPhens = 0;
      size_t commonPhens = 0;
      size_t singlePhens = 0;
      const size_t genomes = tree. root->getLeavesSize ();
      FOR_START (size_t, i, tree. genes, tree. features. size ())
      {
        const Feature& f = tree. features [i];
        f. qc ();
        if (f. genomes == 0)
          optionalPhens++;
        else if (f. genomes == genomes)
        {
          ASSERT (f. gains == 1);
          ASSERT (f. losses == 0);
          commonPhens++;
        }
        else if (f. genomes == 1)
          singlePhens++;
        else   // Non-trivial phenotype
        {
          if (! f. gains)
          {
            f. print (cout);
            ERROR;
          }
          phenGains += f. gains - 1;
          phenLosses += f. losses;
          if (f. gains == 1)
            monoPhens++;
          f. print (cout);
        }
      }
      cout << "# Paraphyletic phenotypes: " << monoPhens   << endl;  // Better: more
      cout << "# Non-paraphyletic gains:  " << phenGains   << endl;  // Better: less
      cout << "# Losses:  "                 << phenLosses  << endl;  // Better: less
      cout << "# Common phenotypes:       " << commonPhens << endl;
      cout << "# Single phenotypes:       " << singlePhens << endl;
      cout << "# Optional phenotypes:     " << optionalPhens << endl;
    }

    
    if (! phen_gains. empty ())
    {
      OFStream of (phen_gains);
     	for (const DiGraph::Node* node : tree. nodes)
     	{
     		const Phyl* phyl = static_cast <const Phyl*> (node);
        FOR_START (size_t, i, tree. genes, tree. features. size ())
          if ((! phyl->feature2parentCore (i) || phyl == tree. root) && phyl->core [i])  
          {
            phyl->getTipName (). saveText (of);
            of << '\t' << tree. features [i]. name << endl;
          }
      }
    }

    
    if (! newick. empty ())
    {
      OFStream os (newick); 
    	os << fixed << setprecision (6);  // PAR
      tree. printNewick (os, true, min_newick_name);
    }
    

    if (! arc_length_stat. empty ())
    {
      OFStream f (arc_length_stat);
      ONumber on (f, 6, false);  // PAR
      tree. printArcLengths (f);
    }


    if (! patrDistFName. empty ())
    {
      const Vector<Tree::Patristic> patrs (tree. getLeafDistances ());
      OFStream f (patrDistFName);
      ONumber on (f, 6, false);  // PAR
      for (const auto patr : patrs)
        f         << patr. leaf1->getName ()
          << '\t' << patr. leaf2->getName ()
          << '\t' << patr. distance 
          << endl;
    }


    if (verbose ())
    {
      tree. setLenGlobal (); 
      tree. setCore ();
      tree. qc ();
      if (verbose ())
      	tree. print (cout);
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


