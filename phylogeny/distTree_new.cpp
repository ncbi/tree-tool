// distTree_new.cpp

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
		: Application ("Find location of new objects in a distance tree")
		{
		  // Input
		  addPositional ("data", "Directory with data");
		  addFlag ("init", "Initialize search");
		  
  	  addKey ("variance", "Dissimilarity variance function: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
  	  addKey ("variance_power", "Power for -variance pow; >= 0", "NaN");
  	  addFlag ("variance_dissim", "Variance is computed off dissimilarities");
  	  addKey ("variance_min", "Min. dissimilarity variance; to be added to the computed variance. If > 0 then all objects are discernible", "0");
  	    	  
		  addKey ("name", "Name of the object");
		  addKey ("dissim", "File of the format: <obj1> <obj2> <dissimilarity>");
		  
		  // Output
		  addKey ("request", "Output file of the format: <obj1> <obj2>");
		  addKey ("leaf", "Output file of the format: <obj_new> <obj1>-<obj2> <leaf_len> <arc_len>");
		}
	
	
	
	void body () const final
  {
	  const string dataDir       = getArg ("data");
	  const bool   init          = getFlag ("init");

	               varianceType    = str2varianceType (getArg ("variance"));  // Global
	               variancePower   = str2real (getArg ("variance_power"));    // Global
	               variance_min    = str2real (getArg ("variance_min"));      // Global
	//const bool   variance_dissim =           getFlag ("variance_dissim");  

	  const string name          = getArg ("name");
	  const string dissimFName   = getArg ("dissim");
	  const string requestFName  = getArg ("request");
	  const string leafFName     = getArg ("leaf");
	   
    if (! isRight (dataDir, "/"))
      throw runtime_error (strQuote (dataDir) + " must end with '/'");

		if (! isNan (variancePower) && varianceType != varianceType_pow)
		  throw runtime_error ("-variance_power requires -variance pow");
		if (isNan (variancePower) && varianceType == varianceType_pow)
		  throw runtime_error ("-variance_power is needed by -variance pow");
		if (variancePower <= 0.0)
		  throw runtime_error ("-variance_power must be positive");
      
    ASSERT (name. empty () == dissimFName.  empty ());
    ASSERT (name. empty () == requestFName. empty ());
    ASSERT (name. empty () == leafFName.    empty ());


    if (verbose ())
    {
      DistTree::printParam (cout);
      cout << endl;
    }

    DistTree tree (dataDir, string (), false, false, false);
    tree. qc ();     

    if (verbose ())
    {
      tree. printInput (cout);
      cout << endl;
    }
    
    if (name. empty ())
    {
      const string newDir (dataDir + "search/");
      FileItemGenerator fig (1, true, newDir);  // PAR
  	  string item;
  	  while (fig. next (item))
      {
        const NewLeaf nl (tree, newDir, item, init);
        nl. qc ();
      }
    }
    else
    {
      const NewLeaf nl (tree, name, dissimFName, leafFName, requestFName, init);
      nl. qc ();
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


