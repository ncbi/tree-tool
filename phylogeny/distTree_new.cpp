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
		}
	
	
	
	void body () const final
  {
	  const string dataDir      = getArg ("data");
	  const bool init           = getFlag ("init");
	               varianceType = str2varianceType (getArg ("variance"));  // Global    
    if (! isRight (dataDir, "/"))
      throw runtime_error (strQuote (dataDir) + " must end with '/'");


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
    
    const string newDir (dataDir + "search/");
    FileItemGenerator fig (1, true, newDir);  // PAR
	  string item;
	  while (fig. next (item))
    {
      const NewLeaf nl (tree, newDir, item, init);
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


