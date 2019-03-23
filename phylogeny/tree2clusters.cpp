// tree2clusters.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{


const string distName = "dist"; 



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Print indiscernibility clusters of a distance tree")
	{
	  addPositional ("input_tree", "File with the tree and arc lengths");
	}



	void body () const final
  {
		const string input_tree = getArg ("input_tree");
				
    
    const DistTree tree (input_tree, string (), string ()); 
    tree. qc (); 
    
    for (const auto& it : tree. name2leaf)
    {
      const Leaf* leaf = it. second;
      if (! leaf->discernible)
        cout << leaf->name 
             << '\t' << static_cast <const DTNode*> (leaf->getParent ()) -> getFirstDecendant () -> getName ()
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


