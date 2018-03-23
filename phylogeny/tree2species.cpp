// tree2species.cpp

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
		: Application ("Find species in a distance tree. Print species under species")
		{
		  addPositional ("input_tree", "Tree file");
		  addPositional ("species_dist", "Max. distance between objects of the same species");
		}



	void body () const final
  {
	  const string input_tree = getArg ("input_tree");
	  const Real species_dist = str2real (getArg ("species_dist"));
    ASSERT (! input_tree. empty ());

    
    DistTree tree (input_tree, string (), string (), false);
    tree. qc ();    
    
    tree. findSpecies (species_dist);
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


