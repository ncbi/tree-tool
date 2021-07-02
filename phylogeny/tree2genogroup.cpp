// tree2genogroup.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Find genogroups of a distance tree
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;
#include "../version.inc"



namespace 
{
  
  
  
struct Genogroup
{
  VectorPtr<Tree::TreeNode> leaves;
    // Distinct
  const Steiner* lca {nullptr};
    // Root of the genogroup subgraph
    // nullptr <=> singleton
  const Leaf* leader {nullptr};
    // !nullptr
    // Prefer complete, well annotated objects !??
  
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
		: Application ("Find genogroups in a distance tree by single linkage clustering")
		{
		  version = VERSION;
		  // Input
		  addPositional ("input_tree", "Tree file");
		  addPositional ("genogroup_dist", "Max. distance between objects of the same genogroup");
		//addKey ("outlier", "Criterion threshold to remove outliers with");  // ??!
		  // Output
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
	 		const DTNode* node = static_cast <const DTNode*> (node_);
	 		if (const Leaf* leaf = node->asLeaf ())
 		    genogroups [node_->getDisjointCluster ()]. leaves << leaf;
	 	}
    for (auto& it : genogroups)
      it. second. finish ();
      
    /* Proof:
         Let \beta be the barrier (= genogroup_dist).
         If d_{ab} \le \beta and d_{cd} \le \beta then
            one of d_{ac}, d_{ad}, d_{bc} and d_{bd} \le \beta, otherwise
            4 \beta < d_{ac} + d_{ad} + d_{bc} + d_{bd} = 2 d_{ab} + 2 d_{cd} \le 4 \beta #.
         => interior nodes of different genogroups do not intersect.
         A genogroups with its interior nodes is a connected subgraph. It has a unique root.
    */
    VectorPtr<Tree::TreeNode> lcas;  lcas. reserve (genogroups. size ());
        // !nullptr
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


