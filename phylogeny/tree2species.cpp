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
		  addKey ("species_table", "File to print: <object> <species>, where <species> is a species representative object");
		}



	void body () const final
  {
	  const string input_tree    = getArg ("input_tree");
	  const Real species_dist    = str2real (getArg ("species_dist"));
	  const string species_table = getArg ("species_table");
    ASSERT (! input_tree. empty ());

    
    DistTree tree (input_tree, string (), string (), string (), false);
    tree. qc ();    
      
    tree. findSpecies (species_dist);
  
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
	  
    // Print species under species
	 	for (DiGraph::Node* node_ : tree. nodes)
	 	{
	 		Tree::TreeNode* node = const_static_cast <Tree::TreeNode*> (node_);
	 		if (Tree::TreeNode* parent = const_cast <Tree::TreeNode*> (node->getParent ()))
	 			if (   clusters [parent->getDisjointCluster ()] > 1  // <=> cluster contains Leaf's
	 				  &&    node   -> getDisjointCluster () 
	 				     != parent -> getDisjointCluster ()
	 				 )
	 				cout << node->getLcaName () << " - " << parent->getLcaName () << endl;
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


