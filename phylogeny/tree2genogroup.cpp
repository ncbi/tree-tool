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
		  addKey ("genogroup_table", "File with lines: <object> <genogroup>, where <genogroup> is a genogroup representative object");
		  addKey ("genogroups", "File with the names of the interior nodes which are genogroup roots");
		  addKey ("genogroup_under_genogroup", "File with lines: <node1 LCA name> <node2 LCA name>, where nodes belong to different genogroup, but node1 is a child of node2");
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

    
    const bool clustersUsed = ! genogroupsFName. empty () || ! genogroup_under_genogroup. empty ();


    DistTree tree (input_tree, string (), string (), false);
    tree. qc ();    
      
    tree. findGenogroups (genogroup_dist);
  
    const size_t nodesSize = tree. nodes. size ();
    unordered_map <const DiGraph::Node*, string> names;  names. rehash (nodesSize);
    unordered_map <const DiGraph::Node*, VectorPtr<Tree::TreeNode>> clusters;  
    if (clustersUsed)
      clusters. rehash (nodesSize);
	 	for (DiGraph::Node* node_ : tree. nodes)
	 	{
	 		const DiGraph::Node* cluster = node_->getDisjointCluster ();
	 		const DTNode* node = static_cast <const DTNode*> (node_);
	 		if (const Leaf* leaf = node->asLeaf ())
	 		{
	 		  if (clustersUsed)
	 		    clusters [cluster] << leaf;
	 			names [cluster] = leaf->name;
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
	  

    if (clustersUsed)
    { 
      VectorPtr<Tree::TreeNode> lcas;  lcas. reserve (nodesSize);
      {
        Tree::LcaBuffer buf;
        for (const auto& it : clusters)
        {
          const VectorPtr<Tree::TreeNode>& cluster = it. second;
          ASSERT (! cluster. empty ());
          if (cluster. size () == 1)
            continue;
          const Tree::TreeNode* lca = Tree::getLca (cluster, buf);
          lcas << lca;
        }
        lcas. sort ();        
        ASSERT (lcas. isUniq ());
      }
      
      if (! genogroupsFName. empty ())
      {
  	  	OFStream f (genogroupsFName);
  	  	for (const Tree::TreeNode* lca : lcas)
  	      f << static_cast <const DTNode*> (lca) -> getLcaName () << endl;
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
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


