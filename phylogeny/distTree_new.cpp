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
		  addPositional ("data", "Directory with data");
		  addFlag ("init", "Initialize search");
		  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
		  addKey ("name", "Name of the object");
		  addKey ("dissim", "File of the format: <obj1> <obj2> <dissimilarity>");
		  addKey ("request", "Output file of the format: <obj1> <obj2>");
		  addKey ("leaf", "Output file of the format: <obj_new> <obj1>-<obj2> <leaf_len> <arc_len>");
		}
	
	
	
	void body () const final
  {
	  const string dataDir      = getArg ("data");
	  const bool init           = getFlag ("init");
	               varianceType = str2varianceType (getArg ("variance"));  // Global   
	  const string name         = getArg ("name");
	  const string dissimFName  = getArg ("dissim");
	  const string requestFName = getArg ("request");
	  const string leafFName    = getArg ("leaf");
	   
    if (! isRight (dataDir, "/"))
      throw runtime_error (strQuote (dataDir) + " must end with '/'");
      
    ASSERT (name. empty () == dissimFName.  empty ());
    ASSERT (name. empty () == requestFName. empty ());
    ASSERT (name. empty () == leafFName.    empty ());


    if (verbose ())
    {
      DistTree::printParam (cout);
      cout << endl;
    }

    DistTree tree (dataDir, false, false, false);
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


