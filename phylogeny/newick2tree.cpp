// newick2tree.cpp

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
	: Application ("Print a tree")
	{
	  // Input
	  addPositional ("input_tree", "Tree in Newick format");
	}



	void body () const final
  {
		const string input_tree = getArg ("input_tree");    

    const DistTree tree (input_tree);
    tree. qc ();         
    tree. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


