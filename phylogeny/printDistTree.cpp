// printDistTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Print a tree made by makeDistTree")
	{
    // Input
	  addPositional ("input_tree", "Tree file");
	  addKey ("data", dmSuff + "-file without " + strQuote (dmSuff) + " to read object comments");
	  addKey ("dissim", "Dissimilarity attribute name in the <data> file");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]); 
	  addKey ("name_match", "File with lines: <name_old> <tab> <name_new>, to replace leaf names");
    // Output
	  addKey ("format", "newick|itree (makeDistTree output)|ASNT (textual ASN.1)", "newick");
	  addFlag ("min_name", "Minimal leaf names for newick");
	  addFlag ("order", "Order subtrees by the number of leaves descending");
	}



	void body () const final
  {
	  const string input_tree     = getArg ("input_tree");
	  const string dataFName      = getArg ("data");
	  const string dissimAttrName = getArg ("dissim");
	               varianceType   = str2varianceType (getArg ("variance"));  // Global    
	  const string name_match     = getArg ("name_match");
	  const string format         = getArg ("format");
  	const bool min_name         = getFlag ("min_name");
  	const bool order            = getFlag ("order");
    ASSERT (! input_tree. empty ());
    if (dataFName. empty () != dissimAttrName. empty ())
      throw runtime_error ("The both data file and the dissimilarity attribute must be present or absent");
    

    DistTree tree (input_tree, dataFName, dissimAttrName, false);
    if (order)
      tree. sort ();
    if (! dataFName. empty ())
      tree. setNodeAbsCriterion ();
    tree. qc ();    
    
    if (! name_match. empty ())
    {
      LineInput f (name_match, 10 * 1024, 1000);  // PAR
      string name_old, name_new;
      while (f. nextLine ())
      {
        name_new = f. line;
        name_old = findSplit (name_new, '\t');
        ASSERT (! name_old. empty ());
        ASSERT (! name_new. empty ());
        if (const Leaf* leaf = findPtr (tree. name2leaf, name_old))
          var_cast (leaf) -> name = name_new;
        else
          throw runtime_error ("Object '" + name_old + "' does not exist");
      }
    }

   	cout << fixed << setprecision (6);  // PAR
    if (format == "newick")
      tree. printNewick (cout, false, min_name);
    else if (format == "itree")
      tree. saveText (cout);
    else if (format == "ASNT")
      tree. printAsn (cout);
    else
      throw runtime_error ("Unknown format " + format);
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


