// replaceDistTree_match.cpp

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
*   Replace a subtree in a "whole" distance tree by matching the subtrees of a "part" tree with the subtrees of the "whole" tree
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{


struct ThisApplication : Application
{
	ThisApplication ()
		: Application ("Replace a subtree in a \"whole\" distance tree by matching the subtrees of a \"part\" tree with the subtrees of the \"whole\" tree")
		{
		  version = VERSION;
		  
		  // Input
		  addPositional ("whole", "Rooted tree where there is a subtree to be replaced by <part>");
		  addPositional ("part", "Rooted tree to replace a subtree in <whole>");
		  
		  // Output
		  addPositional ("out", "Output tree");
		}
	
	
	
	void body () const final
  {
	  const string wholeFName = getArg ("whole");
	  const string partFName  = getArg ("part");
	  const string outFName   = getArg ("out");

	   
    DistTree whole (wholeFName);
    whole. node2deformationPair. clear ();
    whole. qc ();     
    cout << "# Objects in whole tree: " << whole. name2leaf. size () << endl;

    DistTree part (partFName);
    part. setLeaves ();
    part. qc ();     
    cout << "# Objects in part tree: " << part. name2leaf. size () << endl;
      
    Set<const Tree::TreeNode*> whole_leaves; 
    for (const auto& it : whole. name2leaf)
      if (contains (part. name2leaf, it. first))
        whole_leaves << it. second;
    ASSERT (whole_leaves. size () <= part. name2leaf. size ());
    if (whole_leaves. size () < part. name2leaf. size ())
      throw runtime_error ("objects of the part tree are not a subset of the objects of the whole tree");  
    
    const DTNode* root = nullptr;
    {
      const VectorPtr<Tree::TreeNode> roots (whole. leaves2lcas (whole_leaves));
      cout << "Roots:" << endl;
      for (const Tree::TreeNode* node : roots)
        cout << node->getLcaName () << ": " << node->leaves << " objects" << endl;        
      if (roots. size () != 1)
        throw runtime_error ("Multiple roots of the part tree in the whole tree");  
      root = static_cast<const DTNode*> (roots. front ());
    }
    ASSERT (whole. root->leaves == whole. name2leaf. size ());        
    ASSERT (root);
    ASSERT (& root->getDistTree () == & whole);
    ASSERT (root->asSteiner ());
        
    var_cast (root->asSteiner ()) -> replaceSubtree (part);        
    whole. setName2leaf ();
  //whole. setDiscernibles ();  // Requires: optimizable()
    whole. qc ();
        
    {
      OFStream f (outFName);
      whole. saveText (f);
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


