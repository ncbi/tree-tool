// replaceDistTree.cpp

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
*   Add new objects to a distance tree
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
		: Application ("Replace a subtree in a distance tree, the subtree is rerooted")
		{
		  version = VERSION;
		  
		  // Input
		  addPositional ("whole", "Tree where there is a subtree to be replaced by <part>");
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
      

    if (verbose ())
    {
      DistTree::printParam (cout);
      cout << endl;
    }
  

    DistTree whole (wholeFName, noString, noString, noString);
    whole. qc ();     
    cout << "# Objects in whole tree: " << whole. name2leaf. size () << endl;

    DistTree part (partFName, noString, noString, noString);
    part. qc ();     
    cout << "# Objects in part tree: " << part. name2leaf. size () << endl;
      
    Set<const Leaf*> common;  // Leaf's of whole
    for (const auto& it : whole. name2leaf)
      if (contains (part. name2leaf, it. first))
        common << it. second;
    cout << "# Common objects: " << common. size () << endl;
    if (common. empty ())
      throw runtime_error ("No common objects");  // Output whole as the resulk ??
    
    // DTNode::leaves
    var_cast (whole. root) -> setLeaves ();
    ASSERT (whole. root->leaves == whole. name2leaf. size ());

    Vector<Set<const DTNode*>> leaves2nodes;  
      // Index = DTNode::leaves - 1
    leaves2nodes. resize (whole. name2leaf. size ());
    for (const DiGraph::Node* node : whole. nodes)
    {
      const DTNode* dtNode = static_cast <const DTNode*> (node);
      ASSERT (dtNode);
      ASSERT (dtNode->leaves);
      leaves2nodes [dtNode->leaves - 1] << dtNode;
    }    
    
    VectorPtr<DTNode> roots;
    {
      // Bottom-up
      for (Set<const DTNode*>& nodes : leaves2nodes)
      {
        VectorPtr<DTNode> bads; bads. reserve (nodes. size ());
        for (const DTNode* parent : nodes)
        {
          ASSERT (parent);
          bool good = true;
          const VectorPtr<DiGraph::Node> children (parent->getChildren ());
            // Current roots of subtrees
          if (children. empty ())
          {
            ASSERT (parent->leaves == 1);
            ASSERT (parent->asLeaf ());
            if (! common. contains (parent->asLeaf ()))
              good = false;
          }
          else
            for (const DiGraph::Node* child_ : children)
            {
              const DTNode* child = static_cast <const DTNode*> (child_);
              ASSERT (child);
              ASSERT (child->leaves < parent->leaves);
              ASSERT (child->leaves);
              if (! leaves2nodes [child->leaves - 1]. contains (child))
              {
                good = false;
                break;
              }
            }
          if (good)
            for (const DiGraph::Node* child_ : children)
            {
              const DTNode* child = static_cast <const DTNode*> (child_);
              EXEC_ASSERT (leaves2nodes [child->leaves - 1]. erase (child) == 1);
            }
          else
            bads << parent;
        }
        for (const DTNode* bad : bads)
          EXEC_ASSERT (nodes. erase (bad) == 1);
      }
      for (const Set<const DTNode*>& nodes : leaves2nodes)
        for (const DTNode* node : nodes)
          roots << node;
    }      
    ASSERT (! roots. empty ());
    cout << "Roots:" << endl;
    for (const DTNode* node : roots)
      cout << node->getLcaName () << ": " << node->leaves << " objects" << endl;        
    if (roots. size () != 1)
      throw runtime_error ("Multiple roots of the part tree in the whole tree");  // Process each root separately ??
      
    const DTNode* root = roots. front ();
    ASSERT (root);
    ASSERT (& root->getDistTree () == & whole);
    ASSERT (root->asSteiner ());
    
    const double wholeLen =       root->getSubtreeLength ();
    const double partLen  = part. root->getSubtreeLength ();
    ASSERT (wholeLen > 0.0);
    ASSERT (partLen > 0.0);
    const double ratio = wholeLen / partLen;
    QC_ASSERT (ratio > 0.0);
    cout << "Length ratio: " << ratio << endl;
  
  
    // part
    Vector<NewLeaf::Leaf2dissim> leaf2dissims;  leaf2dissims. reserve (common. size ());
    {  
      Tree::LcaBuffer buf;
      for (const Leaf* leaf : common)
      {
        ASSERT (leaf);
        const Tree::TreeNode* lca = nullptr;
        VectorPtr<Tree::TreeNode>& path = Tree::getPath (root, leaf, root, lca, buf);
        ASSERT (lca == root);
        const Real dissim = DistTree::path2prediction (path);
        ASSERT (dissim >= 0.0);
        ASSERT (DM_sp::finite (dissim));
        leaf2dissims << NewLeaf::Leaf2dissim (findPtr (part. name2leaf, leaf->getName ()), dissim / ratio, NaN);
      }
    }    
    ASSERT (leaf2dissims. size () == common. size ());
    
    const NewLeaf nl (part, std::move (leaf2dissims));
    ASSERT (leaf2dissims. empty ());
    nl. qc ();
    cout << nl. location << endl;
    
    part. reroot (var_cast (nl. location. anchor), nl. location. arcLen);
      // nl.location.leafLen is not used


    // whole    
    whole. node2deformationPair. clear ();
    ASSERT (part. root);
    static_cast <const DTNode*> (part. root) -> copySubtree (* var_cast (root->asSteiner ()), ratio);
    
    for (const Leaf* leaf : common)
      whole. removeLeaf (var_cast (leaf), false);

    whole. setName2leaf ();
  //whole. setDiscernibles ();  // Requires: optimizable()
    
    
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


