// tree2genogroup.cpp

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
		: Application ("Find genogroups in a distance tree")
		{
		  addPositional ("input_tree", "Tree file");
		  addPositional ("genogroup_dist", "Max. distance between objects of the same genogroup");
		  addKey ("genogroup_table", "File to print: <object> <genogroup>, where <genogroup> is a genogroup representative object");
		  addKey ("genogroup_under_genogroup", "File to print: <node1 LCA name> <node2 LCA name>, where nodes belong to different genogroup, but node1 is a child of node2");
		}



	void body () const final
  {
	  const string input_tree                = getArg ("input_tree");
	  const Real genogroup_dist              = str2real (getArg ("genogroup_dist"));
	  const string genogroup_table           = getArg ("genogroup_table");
	  const string genogroup_under_genogroup = getArg ("genogroup_under_genogroup");
    ASSERT (! input_tree. empty ());
    ASSERT (genogroup_dist > 0);

    
    DistTree tree (input_tree, string (), string (), false);
    tree. qc ();    
      
    tree. findGenogroups (genogroup_dist);
  
    map <const DiGraph::Node*, VectorPtr<Tree::TreeNode>> clusters;
    map <const DiGraph::Node*, string> names;
	 	for (DiGraph::Node* node_ : tree. nodes)
	 	{
	 		const DiGraph::Node* cluster = node_->getDisjointCluster ();
	 		const DTNode* node = static_cast <const DTNode*> (node_);
	 		if (const Leaf* leaf = node->asLeaf ())
	 		{
	 		  if (! genogroup_under_genogroup. empty ())
	 		    clusters [cluster] << leaf;
	 			names [cluster] = leaf->name;
	 		}
	 	}
	  

    if (! genogroup_under_genogroup. empty ())
    { 
      Set<const Tree::TreeNode*> lcas;
      Tree::LcaBuffer buf;
      size_t n = 0;
      for (const auto& it : clusters)
      {
        const VectorPtr<Tree::TreeNode>& cluster = it. second;
        ASSERT (! cluster. empty ());
        if (cluster. size () == 1)
          continue;
        const Tree::TreeNode* lca = Tree::getLca (cluster, buf);
        lcas << lca;
        n++;
      }
      ASSERT (n == lcas. size ());
      
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

	  
	  if (! genogroup_table. empty ())
	  {
	  	OFStream f (genogroup_table);
		 	for (DiGraph::Node* node_ : tree. nodes)
		 	{
		 		const DTNode* node = static_cast <const DTNode*> (node_);
		 		if (const Leaf* leaf = node->asLeaf ())
		 			f << leaf->name << '\t' << names [node_->getDisjointCluster ()] << endl;
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


