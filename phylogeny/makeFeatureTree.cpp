// makeFeatureTree.cpp

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
*   Make a feature tree
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "featureTree.hpp"
using namespace FeatureTree_sp;
#include "version.inc"



namespace 
{



struct ThisApplication : Application
{
	ThisApplication ()
  	: Application ("Optimize a feature tree")
  	{
  	  version = VERSION;
  	  
  		// Input
  	  addKey ("input_tree", "Input file with the tree");
  	  addKey ("features", "Input directory with features for each genome. Line format: " + Genome::featureLineFormat ());
  	  addFlag ("large", "Featrue files are grouped into subdirectories which are their hash-names (hash<string> % 1000)");
  	  addKey ("input_core", "Input file with root core feature ids");
  	  addFlag ("nominal_singleton_is_optional", "Nominal singleton value means that all values of this nominal attribute are optional for the genome");
  	  addFlag ("prefer_gain", "Prefer gain over loss in maximum parsimony method");
  	    
  	  // Process
  	  addFlag ("use_time", "Use time for MLE, otherwise maximum parsimony method");
  	  addKey ("optim_iter_max", "# Iterations for tree optimization; -1: optimize time only", "0");
  	  addFlag ("save_mem", "Save RAM memory by processing features one by one in the tree. This makes processing slower and restricts functionality");
  	  addKey ("output_core", "Find root, set root core, sort tree and save file with root core feature ids");  	    
  
      // Output	    
  	  addKey ("output_tree", "Output file with the tree");
  	  addKey ("newick", "Output file with the tree in the Newick format");
  	  addFlag ("min_newick_name", "Minimal leaf names in Newick");
  	  addKey ("qual", "Print the summary gain/loss statistics measured by feature consistency, save gain/loss statistcis per feature in the indicated file: +<gaines> -<losses> / <genomes>");
  	  addFlag ("qual_nonredundant", "Non-redundify features by removing identical ones for the summary gain/loss statistics");
  	  addKey ("gain_nodes", "File name to save nodes where features are gained; within a node features are ordered alphabetically");
  	  addKey ("disagreement_nodes", "File name to save nodes with features not gained monophyletically");
  	  addKey ("arc_length_stat", "File with arc length statistics in format " + Tree::printArcLengthsColumns ());
  	  addKey ("patr_dist", "File with patristic distances in format: <leaf name1> <leaf name2> <distance>, where <leaf name1> < <leaf name2>");
  	}



	void body () const final
  {
		const string input_tree         = getArg ("input_tree");
		const string feature_dir        = getArg ("features");
		const bool   large              = getFlag ("large");
		const string input_core         = getArg ("input_core");
		const bool   nominal_singleton_is_optional = getFlag ("nominal_singleton_is_optional");
		const bool   prefer_gain        = getFlag ("prefer_gain");

		const bool   use_time           = getFlag ("use_time");  
		const int    optim_iter_max     = str2<int> (getArg ("optim_iter_max"));
		const bool   save_mem           = getFlag ("save_mem");
		const string output_core        = getArg ("output_core");

		const string output_tree        = getArg ("output_tree");
		const string newick             = getArg ("newick");
		const bool   min_newick_name    = getFlag ("min_newick_name");
		const string qual               = getArg ("qual");  
		const bool   qual_nonredundant  = getFlag ("qual_nonredundant");
		const string gain_nodes         = getArg ("gain_nodes");  
		const string disagreement_nodes = getArg ("disagreement_nodes");  
		const string arc_length_stat    = getArg ("arc_length_stat");
		const string patrDistFName      = getArg ("patr_dist");

		if (save_mem && (optim_iter_max || use_time))
		  throw runtime_error ("-save_mem cannot be used with -optim_iter_max or -use_time");
		if (! input_core. empty () && ! output_core. empty () && input_core == output_core)
		  throw runtime_error ("-input_core file cannot be the same as -output_core file");
		if (optim_iter_max < -1)
		  throw runtime_error ("-optim_iter_max cannot be less than -1");
		if (qual_nonredundant && qual. empty ())
		  throw runtime_error ("if -qual_nonredundant then -qual cannot be empty");
		
		
    const Chronometer_OnePass cop ("Total");  
		
		
    FeatureTree tree (input_tree, feature_dir, large, input_core, nominal_singleton_is_optional, prefer_gain, save_mem);
    tree. printInput (cout);
    tree. qc ();    
    

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
          tree. setCore ();
    	    tree. saveFile (output_tree);  
          if (iter >= optim_iter_max)
          	break;
          iter++;
    	    cout << endl;
    	    cout << "Iter = " << iter << endl;
          if (! tree. optimize ())
          	break;	    
        }
        tree. setCore ();
        tree. saveFile (output_tree);  
      }
      else
      {
   	    cout << endl;
  	  	cout << "Optimizing time ..." << endl;
        tree. optimizeTime ();
      }
    }
    
    
    if (! output_core. empty ())
    {
 	    cout << endl;
      if (tree. allTimeZero)
      {
        // 63 sec./50K genomes
        cout << "Old root: " << tree. root->getLcaName () << endl;
        string newRootName;
        if (save_mem)
        {
          if (const Species* newRoot = tree. findRoot ())
            newRootName = newRoot->getLcaName ();
        }
        else
        {
         	for (const DiGraph::Node* node : tree. nodes)
         		if (const Species* sp = static_cast <const Phyl*> (node) -> asSpecies ())
         		  var_cast (sp) -> setMiddleCoreSize ();
         	newRootName = tree. changeRoot ();
        }
        cout << "New root: " << newRootName << endl;
        tree. sort ();
      }
      else
      {
        ASSERT (! save_mem);
        cout << "Adjusting root core ..." << endl;
        tree. setCore ();
        size_t coreChange [2/*core2nonCore*/];
        tree. resetSuperRootCore (coreChange);
        cout << "# Core to non-core: " << coreChange [true]  << endl;
        cout << "# Non-core to core: " << coreChange [false] << endl;
      }
      // Output
      if (save_mem)
        cout << "FILE WITH ROOT CORE FEATURE IDS IS NOT SAVED!" << endl;
      else
        tree. saveSuperRootCore (output_core);
    }


    if (! save_mem)
      tree. setStats ();
    tree. qc ();


    // Output
    tree. saveFile (output_tree);


    if (! qual. empty ())
    {    
    	// Independent of tree.root

    	Vector<Feature> features = tree. features;
	    if (qual_nonredundant)
	    {
	    	features. sort (Feature::statLess);
	    	features. uniq (Feature::statEqual);
	    }
	
      OFStream out (qual);
      // Input: Feature::{genomes,gains,losses}
      cout << endl;
      cout << "Feature statistics:" << endl;
      size_t monophyletics = 0;     // In an unrooted tree
      size_t nonMonophyletics = 0;  // In an unrooted tree
      size_t extraMutations = 0;
      size_t optionals = 0;
      size_t commons = 0;  // Feature may be optional
      size_t singles = 0;  // Feature may be optional
      const size_t genomes = tree. root->getLeavesSize ();
      for (const Feature& f : features)
      {
        f. qc ();
        if (f. genomes == 0)
          optionals++;
        else if (f. genomes == genomes)
        {
          ASSERT (f. realGains () == 1);
          ASSERT (f. losses. empty ());
          commons++;
        }
        else if (f. genomes == 1)
          singles++;
        else   // Non-trivial features
        {
        	ASSERT (f. mutations () >= 1);
          if (f. mutations () == 1)
            monophyletics++;
          else
          {
            extraMutations += f. mutations () - 1;
            nonMonophyletics++;
          }
          f. print (out);
        }
      }
      const ONumber on (cout, 2, false);  // PAR
      cout << "# Monophyletic features:          " << monophyletics    << " (" << (Real) monophyletics    / (Real) genomes << ") ^" << endl;  
      cout << "# Non-monophyletic features:      " << nonMonophyletics << " (" << (Real) nonMonophyletics / (Real) genomes << ") V" << endl;  
      cout << "# Non-monophyletic disagreements: " << extraMutations   << " (" << (Real) extraMutations   / (Real) genomes << ") V !" << endl;  // = Taxonomy miscongruence for taxonomy features
      cout << "# Common features:                " << commons << endl;
      cout << "# Single features:                " << singles << endl;
      cout << "# Optional features:              " << optionals << endl;
      ASSERT (optionals + commons + singles + monophyletics + nonMonophyletics == features. size ());
      if (! qual_nonredundant && ! use_time && (size_t) round (tree. len) - tree. globalSingletonsSize != monophyletics + nonMonophyletics + extraMutations + singles)
      {
      #if 0
        if (! save_mem)
        {
          size_t s = 0;
          FFOR (size_t, i, features. size ())
          {
            size_t feature_paraphyletics = 0;  
            size_t feature_nonParaphyletics = 0;  
            size_t feature_extraMutations = 0;
            size_t feature_optionals = 0;
            size_t feature_commons = 0;  
            size_t feature_singles = 0;  
            const Feature& f = features [i];
            f. qc ();
            if (f. genomes == 0)
              feature_optionals++;
            else if (f. genomes == genomes)
            {
              ASSERT (f. realGains () == 1);
              ASSERT (f. losses. empty ());
              feature_commons++;
            }
            else if (f. genomes == 1)
              feature_singles++;
            else   // Non-trivial features
            {
            	ASSERT (f. mutations () >= 1);
              if (f. mutations () == 1)
                feature_paraphyletics++;
              else
              {
                feature_extraMutations += f. mutations () - 1;
                feature_nonParaphyletics++;
              }
            }
            const size_t feature_len = (size_t) round (tree. feature2treeLength (i));
            if (feature_len != feature_paraphyletics + feature_nonParaphyletics + feature_extraMutations + feature_singles)
            {
              cout << f;
              cout        << feature_len 
                   << ' ' << feature_paraphyletics 
                   << ' ' << feature_nonParaphyletics 
                   << ' ' << feature_extraMutations 
                   << ' ' << feature_singles
                   << ' ' << f. gains. contains (static_cast <const Phyl*> (tree. root))
                   << ' ' << f. rootGain
                   << endl;
              ERROR;
            }
            s += feature_len;
          }
          cout        << tree. len 
               << ' ' << tree. globalSingletonsSize 
               << ' ' << monophyletics 
               << ' ' << nonMonophyletics 
               << ' ' << extraMutations 
               << ' ' << singles
               << ' ' << s
               << endl;
        }
      #endif
        ERROR;
      }
    }

    
    if (! gain_nodes. empty ())
    {
      OFStream f (gain_nodes);
      for (const Feature& feature : tree. features)
      {        
        for (const Phyl* phyl : feature. gains)
          f << phyl->getLcaName () << '\t' << feature. name << endl;          
        if (feature. rootGain)
          f << tree. root->getLcaName () << '\t' << feature. name << endl;          
      }
  	 	for (const DiGraph::Node* node : tree. nodes)
  	 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
  	 		  for (const Feature::Id& featureId : g->singletons)
  	 		    if (! Feature::nominalSingleton (featureId))
  	 		      f << g->getLcaName () << '\t' << featureId << endl;
    }

    
    if (! disagreement_nodes. empty ())
    {
      OFStream f (disagreement_nodes);
      for (const Feature& feature : tree. features)
        if (! feature. monophyletic ())
        {
          for (const Phyl* phyl : feature. gains)
            f << phyl->getLcaName () << '\t' << feature. name << "\tgain" << endl;
          for (const Phyl* phyl : feature. losses)
            f << phyl->getLcaName () << '\t' << feature. name << "\tloss" << endl;
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
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}


