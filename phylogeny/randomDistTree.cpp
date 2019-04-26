// randomDistTree.cpp

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
	: Application ("Print a random distance tree")
	{
	  // Input
	  addPositional ("branch_prob", "Probability to expand a branch");
	  addPositional ("leaf_num_max", "Max. number of leaves");
	}



	void body () const final
  {
		const Real branch_prob    = str2<Prob> (getArg ("branch_prob"));
		const size_t leaf_num_max = str2<size_t> (getArg ("leaf_num_max"));
		ASSERT (isProb (branch_prob));
		ASSERT (branch_prob < 1.0);
		ASSERT (branch_prob > 0.0);
		ASSERT ((bool) leaf_num_max);
    

    DistTree tree (branch_prob, leaf_num_max);
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


