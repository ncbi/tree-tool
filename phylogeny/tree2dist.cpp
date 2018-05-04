// tree2dist.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{


const string distName = "dist"; 



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Print a " + dmSuff + "-file with an attribute '" + distName + "'")
	{
	  // Input
	  addKey ("input_tree", "File with the tree and arc lengths");
	  addKey ("branch_prob", "Probability to expand a branch", "0");
	  addKey ("leaf_num_max", "Max. number of leaves", "0");
	  addKey ("noise", "dost += Normal(0,noise)", "0");
	  addFlag ("sqrt", "dist = sqrt(dist)");
	//seed ??
    // Output
	  addKey ("output_tree", "Output tree file");
	}



	void body () const final
  {
		const string input_tree   = getArg ("input_tree");
		const Real branch_prob    = str2<Prob> (getArg ("branch_prob"));
		const size_t leaf_num_max = str2<size_t> (getArg ("leaf_num_max"));
		const Real noise          = str2<Real> (getArg ("noise"));
		const bool sqrtP          = getFlag ("sqrt");
		const string output_tree  = getArg ("output_tree");
		ASSERT (isProb (branch_prob));
		ASSERT (branch_prob < 1);
		ASSERT (input_tree. empty () == (bool) branch_prob);
		ASSERT ((bool) leaf_num_max == (bool) branch_prob);
		ASSERT (noise >= 0);    
    

    Common_sp::AutoPtr<DistTree> tree;
    if (input_tree. empty ())
      tree = new DistTree (branch_prob, leaf_num_max);
    else
      tree = new DistTree (input_tree, string (), string (), false); 
    ASSERT (tree);
    tree->qc ();     
      
  //tree->printInput (cout);   
    

    Dataset ds;
    
    VectorPtr<Leaf> leaves;
    for (const auto it : tree->name2leaf)
    {
      ds. appendObj (it. first);
      ASSERT (it. second->graph);
      leaves << it. second;
    }
    ASSERT (ds. objs. size () == leaves. size ());
    
    auto distAttr = new PositiveAttr2 (distName, ds);
    Normal norm;
    norm. setParam (0, noise); 
    FOR (size_t, row, ds. objs. size ())
    {
      distAttr->put (row, row, 0);
      FOR (size_t, col, row)
      {
        const Tree::TreeNode* lca = nullptr;
        const VectorPtr<Tree::TreeNode> path (tree->getPath (leaves [row], leaves [col], nullptr, lca));
        Real dist = tree->path2prediction (path);
        if (noise)
          dist += norm. rand ();
        if (sqrtP)
          dist = sqrt (dist);  
        distAttr->putSymm (row, col, dist);
      }
    }

      
    ds. saveText (cout);    
    
  
    if (! output_tree. empty ())
      tree->saveFile (output_tree);  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


