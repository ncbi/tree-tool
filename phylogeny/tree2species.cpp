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
		: Application ("Find species in a distance tree")
		{
		  addPositional ("input_tree", "Tree file");
		  addPositional ("species_dist", "Max. distance between objects of the same species");
		  addKey ("species_table", "File to print: <object> <species>, where <species> is a species representative object");
		  addKey ("species_under_species", "File to print: <node1 LCA name> <node2 LCA name>, where nodes belong to different species, but node1 is a child of node2");
		}



	void body () const final
  {
	  const string input_tree            = getArg ("input_tree");
	  const Real species_dist            = str2real (getArg ("species_dist"));
	  const string species_table         = getArg ("species_table");
	  const string species_under_species = getArg ("species_under_species");
    ASSERT (! input_tree. empty ());
    ASSERT (species_dist > 0);

    
    DistTree tree (input_tree, string (), string (), false);
    tree. qc ();    
      
    tree. findSpecies (species_dist, ! species_under_species. empty ());
  
    map <const DiGraph::Node*, size_t> clusters;
    map <const DiGraph::Node*, string> names;
	 	for (DiGraph::Node* node_ : tree. nodes)
	 	{
	 		const DiGraph::Node* cluster = node_->getDisjointCluster ();
	 		clusters [cluster] ++;
	 		const DTNode* node = static_cast <const DTNode*> (node_);
	 		if (const Leaf* leaf = node->asLeaf ())
	 			names [cluster] = leaf->name;
	 	}
	  
    if (! species_under_species. empty ())
    { 
      // Traverse the whole path from child to root ??!
	  	OFStream f (species_under_species);
  	 	for (DiGraph::Node* node_ : tree. nodes)
  	 	{
  	 		Tree::TreeNode* node = const_static_cast <Tree::TreeNode*> (node_);
  	 		if (Tree::TreeNode* parent = const_cast <Tree::TreeNode*> (node->getParent ()))
  	 			if (   clusters [parent->getDisjointCluster ()] > 1  // <=> cluster contains Leaf's
  	 				  &&    node   -> getDisjointCluster () 
  	 				     != parent -> getDisjointCluster ()
  	 				 )
  	 				f << node->getLcaName () << " - " << parent->getLcaName () << endl;
  	  }
  	}
	  
	  if (! species_table. empty ())
	  {
	  	OFStream f (species_table);
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


