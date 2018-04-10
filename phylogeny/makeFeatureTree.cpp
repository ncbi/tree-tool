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
  	  addKey ("features", "Input directory with features for each genome. Line format: " + Genome::geneLineFormat ());
  	  addKey ("input_core", "Input file with root core feature ids");
  	    
  	  // Optimization
  	  addFlag ("use_time", "Use time for MLE, otherwise parsimony method");
  	  addKey ("optim_iter_max", "# Iterations for tree optimization; -1: optimize time only", "0");
  	  addKey ("output_core", "Find root, set root core and save file with root core feature ids");  	    
  
      // Output	    
  	  addKey ("output_tree", "Output file with the tree");
  	  addKey ("report_feature", "Feature name to report in the output tree");
  	//addFlag ("set_node_ids", "Set node ids in the output tree");  
  	  addKey ("newick", "Output file with the tree in the Newick format");
  	  addFlag ("min_newick_name", "Minimal leaf names in Newick");
  	  addKey ("qual", "Print the summary gain/loss statistics measured by feature consistency, save gain/loss statistcis per feature in the indicated file: +<gaines> -<losses> / <genomes>");
  	  addKey ("gain_nodes", "File name to save nodes where features are gained");
  	  addKey ("disagreement_nodes", "File name to save nodes where features are not gained monophyletically");
  	  addKey ("arc_length_stat", "File with arc length statistics in format " + Tree::printArcLengthsColumns ());
  	  addKey ("patr_dist", "File with patristic distances in format: <leaf name1> <leaf name2> <distance>, where <leaf name1> < <leaf name2>");
  	}



	void body () const final
  {
		const string input_tree         = getArg ("input_tree");
		const string feature_dir        = getArg ("features");
		const string input_core         = getArg ("input_core");

		const bool use_time             = getFlag ("use_time");  
		const int optim_iter_max        = str2<int> (getArg ("optim_iter_max"));

		const string output_tree        = getArg ("output_tree");
		const string report_feature     = getArg ("report_feature");
	//const bool set_node_ids         = getFlag ("set_node_ids");  
		const string output_core        = getArg ("output_core");
		const string newick             = getArg ("newick");
		const bool min_newick_name      = getFlag ("min_newick_name");
		const string qual               = getArg ("qual");  
		const string gain_nodes         = getArg ("gain_nodes");  
		const string disagreement_nodes = getArg ("disagreement_nodes");  
		const string arc_length_stat    = getArg ("arc_length_stat");
		const string patrDistFName      = getArg ("patr_dist");

		IMPLY (! input_core. empty () && ! output_core. empty (), input_core != output_core);
	//ASSERT (save_feature. empty () || save_features. empty ());
	//ASSERT (species_id > 0);
		ASSERT (optim_iter_max >= -1);
		
		
    FeatureTree tree (/*species_id,*/ input_tree, feature_dir, input_core);
    tree. printInput (cout);
    tree. qc ();    
    
    if (! report_feature. empty ())
    {
      tree. reportFeature = tree. findFeature (report_feature);
      if (tree. reportFeature == NO_INDEX)
        throw runtime_error ("Feature " + report_feature + " is not found");
    }
    

   	if (use_time && tree. allTimeZero)
    {
   		tree. useTime (input_core);
      cout << endl;
      tree. printInput (cout);
    }
    
    
    if (optim_iter_max)  
    {
      if (optim_iter_max > 0)  
      {
  	  	cerr << "Optimizing ..." << endl;
        int iter = 0; 
        for (;;)
        {
    	    cout << endl;
    	    tree. dump (output_tree /*, true*/);  // Redundant
          if (iter >= optim_iter_max)
          	break;
          iter++;
    	    cout << "Iter = " << iter << endl;
          if (! tree. optimize ())
          	break;	    
        }
        tree. dump (output_tree /*, true*/);  // Redundant
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
        cout << "Old root: " << tree. root->getLcaName () << endl;
        cout << "New root: " << tree. findRoot () << endl;
      }
      else
      {
        cout << endl;
        cout << "Adjusting root core ..." << endl;
        tree. setCore ();
        size_t coreChange [2/*core2nonCore*/];
        tree. resetSuperRootCore (coreChange);
        cout << "# Core to non-core: " << coreChange [true]  << endl;
        cout << "# Non-core to core: " << coreChange [false] << endl;
      }
      // Output
      tree. saveSuperRootCore (output_core);
    }

    tree. qc ();


    // Output
    tree. dump (output_tree/*, set_node_ids*/);


    if (! qual. empty ())
    {    
      OFStream out (qual);
      // Input: Feature::{genomes,gains,losses}
      cout << endl;
      cout << "Feature statistics:" << endl;
      size_t monos = 0;  // Monophyletic
      size_t gains = 0;
      size_t losses = 0;
      size_t optionals = 0;
      size_t commons = 0;
      size_t singles = 0;
      const size_t genomes = tree. root->getLeavesSize ();
      FFOR (size_t, i, tree. features. size ())
      {
        const Feature& f = tree. features [i];
        f. qc ();
        if (f. genomes == 0)
          optionals++;
        else if (f. genomes == genomes)
        {
          ASSERT (f. gains == 1);
          ASSERT (f. losses == 0);
          commons++;
        }
        else if (f. genomes == 1)
          singles++;
        else   // Non-trivial features
        {
          if (f. gains <= 1)
            monos++;
          else
            gains += f. gains - 1;
          losses += f. losses;
          //
          f. print (out);
        }
      }
      cout << "# Paraphyletic features:  " << monos   << " ^" << endl;  // Better: more
      cout << "# Non-paraphyletic gains: " << gains   << " v" << endl;  // Better: less
      cout << "# Losses:                 " << losses  << " v" << endl;  // Better: less
      cout << "# Common features:        " << commons << endl;
      cout << "# Single features:        " << singles << endl;
      cout << "# Optional features:      " << optionals << endl;
      cout << endl;
      cout << "Feature disagreement: " << gains + losses << endl;  // Better: less  
    }

    
    if (! gain_nodes. empty ())
    {
      OFStream f (gain_nodes);
     	for (const DiGraph::Node* node : tree. nodes)
     	{
     		const Phyl* phyl = static_cast <const Phyl*> (node);
        FFOR (size_t, i, tree. features. size ())
          if ((! phyl->feature2parentCore (i) || phyl == tree. root) && phyl->core [i])  
            f << phyl->getLcaName () << '\t' << tree. features [i]. name << endl;
      }
    }

    
    if (! disagreement_nodes. empty ())
    {
      OFStream f (disagreement_nodes);
     	for (const DiGraph::Node* node : tree. nodes)
     	{
     		const Phyl* phyl = static_cast <const Phyl*> (node);
        FFOR (size_t, i, tree. features. size ())
          if ((! phyl->feature2parentCore (i) || phyl == tree. root) == phyl->core [i]
          	  && ! tree. features [i]. monophyletic ()
          	 )
            f << phyl->getLcaName () << '\t' << tree. features [i]. name << '\t' << phyl->core [i] << endl;
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


