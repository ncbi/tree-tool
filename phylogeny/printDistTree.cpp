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
	  addKey ("data", dmSuff + "-file without \"" + dmSuff + "\"");
	  addKey ("dissim", "Dissimilarity attribute name in the <data> file");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);  
    // Output
	  addKey ("format", "newick|ASN", "newick");
	  addFlag ("min_name", "Minimal leaf names");
	}



	void body () const
  {
	  const string input_tree     = getArg ("input_tree");
	  const string dataFName      = getArg ("data");
	  const string dissimAttrName = getArg ("dissim");
	               varianceType   = str2varianceType (getArg ("variance"));  // Global    	  
	  const string format         = getArg ("format");
  	const bool min_name         = getFlag ("min_name");
    ASSERT (! input_tree. empty ());
    if (dataFName. empty () != dissimAttrName. empty ())
      throw runtime_error ("The both data file and the dissimilarity attribute must be present or absent");
    

    // Reading dissimFName is slow ??!
    // Tree format should contain all quality attributes ??
    DistTree tree (input_tree, dataFName, dissimAttrName, false);
    tree. sort ();
    if (! dataFName. empty ())
      tree. setLeafAbsCriterion ();
    tree. qc ();    

   	cout << fixed << setprecision (6);  // PAR
    if (format == "newick")
      tree. printNewick (cout, false, min_name);
    else if (format == "ASN")
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


