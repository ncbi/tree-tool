// replaceDistTree_reroot.cpp

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
*   Replace a subtree in a \"whole\" distance tree with a rerooted \"part\" tree 
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
		: Application ("Replace a subtree in a \"whole\" distance tree with a rerooted \"part\" tree")
		{
		  version = VERSION;
		  
		  // Input
		  addPositional ("whole", "Rooted tree where there is a subtree to be replaced by <part>");
		  addPositional ("part", "Tree to replace a subtree in <whole>");
		  
  	  addKey ("variance", "Dissimilarity variance function for <part>: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
  	  addKey ("variance_power", "Power for -variance pow for <part>; >= 0", "NaN");
  	  addKey ("variance_min", "Min. dissimilarity variance for <part>", "0");
  	  
		  // Output
		  addPositional ("out", "Output tree");
		}
	
	
	
	void body () const final
  {
	  const string wholeFName    = getArg ("whole");
	  const string partFName     = getArg ("part");
	  const string outFName      = getArg ("out");

	               varianceType  = str2varianceType (getArg ("variance"));  // Global
	               variancePower = str2real (getArg ("variance_power"));    // Global
	               variance_min  = str2real (getArg ("variance_min"));      // Global
	               
	   
		if (! isNan (variancePower) && varianceType != varianceType_pow)
		  throw runtime_error ("-variance_power requires -variance pow");
		if (isNan (variancePower) && varianceType == varianceType_pow)
		  throw runtime_error ("-variance_power is needed by -variance pow");
		if (variancePower <= 0.0)
		  throw runtime_error ("-variance_power must be positive");
      

    DistTree whole (DissimParam (), wholeFName, noString, noString, noString);
    whole. node2deformationPair. clear ();
    whole. qc ();     
    cout << "# Objects in whole tree: " << whole. name2leaf. size () << endl;

    DistTree part (DissimParam (), partFName, noString, noString, noString);
    part. qc ();     
    cout << "# Objects in part tree: " << part. name2leaf. size () << endl;
      
    Set<const Tree::TreeNode*> whole_leaves;  
    for (const auto& it : whole. name2leaf)
      if (contains (part. name2leaf, it. first))
        whole_leaves << it. second;
    cout << "# Common objects: " << whole_leaves. size () << endl;
    if (whole_leaves. empty ())
      throw runtime_error ("No common objects");  // Output whole as the result ??
    
    const DTNode* root = nullptr;
    {
      const VectorPtr<Tree::TreeNode> roots (whole. leaves2lcas (whole_leaves));
      ASSERT (whole. root->leaves == whole. name2leaf. size ());        
      cout << "Roots:" << endl;
      for (const Tree::TreeNode* node : roots)
        cout << node->getLcaName () << ": " << node->leaves << " objects" << endl;        
      if (roots. size () != 1)
        throw runtime_error ("Multiple roots of the part tree in the whole tree");  // Process each root separately ??
      root = static_cast<const DTNode*> (roots. front ());
    }
    ASSERT (root);
    ASSERT (& root->getDistTree () == & whole);
    ASSERT (root->asSteiner ());
    
    // ratio
    const double wholeLen =       root->getSubtreeLength ();
    const double partLen  = part. root->getSubtreeLength ();
    QC_ASSERT (wholeLen > 0.0);
    QC_ASSERT (partLen > 0.0);
    const double ratio = wholeLen / partLen;
    ASSERT (ratio > 0.0);
    cout << "Length ratio: " << ratio << endl;
  
    NewLeaf::Location loc;
    {
      Vector<NewLeaf::Leaf2dissim> leaf2dissims;  leaf2dissims. reserve (whole_leaves. size ());
      {  
        Tree::LcaBuffer buf;
        for (const Tree::TreeNode* leaf_ : whole_leaves)
        {
          const Leaf* leaf = static_cast<const Leaf*> (leaf_);
          ASSERT (leaf);
          const Tree::TreeNode* lca = nullptr;
          const VectorPtr<Tree::TreeNode>& path = Tree::getPath (root, leaf, root, lca, buf);
          ASSERT (lca == root);
          const Real dissim = DistTree::path2prediction (path);
          ASSERT (dissim >= 0.0);
          ASSERT (DM_sp::finite (dissim));
          leaf2dissims << NewLeaf::Leaf2dissim (findPtr (part. name2leaf, leaf->getName ()), dissim / ratio, NaN);
        }
      }    
      ASSERT (leaf2dissims. size () == whole_leaves. size ());      
      const NewLeaf nl (part, std::move (leaf2dissims));
      ASSERT (leaf2dissims. empty ());
      nl. qc ();
      loc = nl. location;
    }      
    loc. qc ();
    cout << loc << endl;
    
    // part
    part. reroot (var_cast (loc. anchor), loc. arcLen);

    // whole    
    ASSERT (part. root);
    const Steiner* st = static_cast<const DTNode*> (part. root) -> asSteiner ();
    ASSERT (st);
    st->copySubtree (* var_cast (root->asSteiner ()), ratio);
    var_cast (root) -> len +=  loc. leafLen * ratio;
    for (const Tree::TreeNode* leaf : whole_leaves)
      whole. removeLeaf (const_static_cast<Leaf*> (leaf), false);
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


