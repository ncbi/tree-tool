// printDistTree.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{



struct ThisApplication : DistTreeApplication
{
	ThisApplication ()
	: DistTreeApplication ("Print a tree made by makeDistTree")
	{
  #if 0
    // Input
	  addPositional ("input_tree", "Tree file");
	  addKey ("dist", dmSuff + "-file without \"" + dmSuff + "\"");
	  addKey ("attr", "Dissimilarity attribute name", "dist");
	  {
  	  string varianceTypes;
  	  for (const string& s : varianceTypeNames)
  	    varianceTypes += " " + s;
  	  addKey ("variance", "Dissimilarity variance: " + varianceTypes, "exp");
  	}
  #endif
    // Output
	  addKey ("format", "newick|ASN", "newick");
	  addFlag ("min_name", "Minimal leaf names");
	}



	void body () const
  {
  #if 1
    const DistTreeParam dtp (*this);
  #else
		const string input_tree     = getArg ("input_tree");
    const string dataFName      = getArg ("data");
	  const string dissimAttrName = getArg ("dissim");
  	             varianceType = str2varianceType (getArg ("variance"));
	#endif
	  
	  const string format       = getArg ("format");
  	const bool min_name       = getFlag ("min_name");
  //ASSERT (! dtp. input_tree. empty ());
    

    // Reading dissimFName is slow ??!
    // Tree format should contain all quality attributes ??
    DistTree tree (dtp. input_tree, dtp. dataFName, dtp. dissimAttrName, false);
    tree. sort ();
    if (! dtp. dataFName. empty ())
      tree. setLeafAbsCriterion ();
    tree. qc ();    

   	cout << fixed << setprecision (6);  // PAR
    if (format == "newick")
      tree. printNewick (cout, false, min_name);
    else if (format == "ASN")
      tree. printAsn (cout);
    else
      ERROR_MSG ("Unknown format " + format);
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


