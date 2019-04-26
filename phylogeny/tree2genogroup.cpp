// tree2genogroup.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{
  
  
  
struct Genogroup
{
  VectorPtr<Tree::TreeNode> leaves;
    // Distinct
  const Steiner* lca {nullptr};
  const Leaf* leader {nullptr};
  
  Genogroup () = default;
    
  void finish ()
    { ASSERT (! leaves. empty ());
      if (leaves. size () == 1)
        leader = static_cast <const DTNode*> (leaves [0]) -> asLeaf ();
      else
      { Tree::LcaBuffer buf;
        const Tree::TreeNode* lca_ = Tree::getLca (leaves, buf);
        EXEC_ASSERT (lca = static_cast <const DTNode*> (lca_) -> asSteiner ());
        EXEC_ASSERT (leader = static_cast <const DTNode*> (lca->getFirstDecendant ()) -> asLeaf ());
      }
      ASSERT (leader);
    }
};



struct ThisApplication : Application
{
	ThisApplication ()
		: Application ("Find genogroups in a distance tree")
		{
		  addPositional ("input_tree", "Tree file");
		  addPositional ("genogroup_dist", "Max. distance between objects of the same genogroup");
		  addKey ("genogroup_table", "File with lines: <object> <genogroup leader>");
		  addKey ("genogroups", "File with the names of the interior nodes which are genogroup roots");
		  addKey ("genogroup_under_genogroup", "File with lines: <node1 LCA name> <node2 LCA name>, where nodes belong to different genogroups, but node1 is a child of node2");
		}



	void body () const final
  {
	  const string input_tree                = getArg ("input_tree");
	  const Real genogroup_dist              = str2real (getArg ("genogroup_dist"));
	  const string genogroup_table           = getArg ("genogroup_table");
	  const string genogroupsFName           = getArg ("genogroups");
	  const string genogroup_under_genogroup = getArg ("genogroup_under_genogroup");
    ASSERT (! input_tree. empty ());
    ASSERT (genogroup_dist > 0);

    
    DistTree tree (input_tree, string (), string (), string());
    tree. sort ();  // For Genogroup::leader
    tree. qc ();    
      
    tree. findGenogroups (genogroup_dist);
  
    unordered_map <const DiGraph::Node* /*genogroup cluster*/, Genogroup> genogroups;  
    genogroups. rehash (tree. nodes. size ());
	 	for (DiGraph::Node* node_ : tree. nodes)
	 	{
	 		const DiGraph::Node* cluster = node_->getDisjointCluster ();
	 		const DTNode* node = static_cast <const DTNode*> (node_);
	 		if (const Leaf* leaf = node->asLeaf ())
 		    genogroups [cluster]. leaves << leaf;
	 	}
    for (auto& it : genogroups)
      it. second. finish ();
      
    VectorPtr<Tree::TreeNode> lcas;  lcas. reserve (genogroups. size ());
    for (const auto& it : genogroups)
      if (const Steiner* lca = it. second. lca)
        lcas << lca;
    lcas. sort ();
    ASSERT (lcas. isUniq ());

	  
	  // Output
	  if (! genogroup_table. empty ())
	  {
	  	OFStream f (genogroup_table);
		 	for (const auto& it : tree. name2leaf)
		 	{
		 		const Leaf* leaf = it. second;
	 		  f << leaf->name << '\t' << genogroups [var_cast (leaf) -> getDisjointCluster ()]. leader->name << endl;
		 	}
	  }
	      
    if (! genogroupsFName. empty ())
    {
	  	OFStream f (genogroupsFName);
      for (const Tree::TreeNode* lca : lcas)
        f << lca->getLcaName () << endl;
	  }

    if (! genogroup_under_genogroup. empty ())
    {
	  	OFStream f (genogroup_under_genogroup);
	  	for (const Tree::TreeNode* lca : lcas)
	  	{
	  	  const Tree::TreeNode* parent = lca->getParent ();
	  	  while (parent)
	  	  {
	  	    if (lcas. contains (parent))
	  	    {
	  	      f << static_cast <const DTNode*> (lca)    -> getLcaName () << " - " 
	  	        << static_cast <const DTNode*> (parent) -> getLcaName () 
	  	        << endl;
	  	      break;
	  	    }
	  	    parent = parent->getParent ();
	  	  }
	  	}
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


