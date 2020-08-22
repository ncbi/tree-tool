// asnt2tree.cpp

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
*   Convert a tree in textual ASN.1 format to internal format
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "distTree.hpp"
#include "asn.hpp"
using namespace DistTree_sp;



namespace 
{


struct TreeAsn : Asn_sp::Asn
{
  DistTree& tree;
  map<size_t/*id*/,Steiner*> id2steiner;
  

  TreeAsn (const string &fName,
           DistTree &tree_arg)
    : Asn (fName)
    , tree (tree_arg)
    {
      ASSERT (tree. nodes. empty ());
    }

    
  void processNode () final
    { 
      constexpr size_t labelNum = 0;
    	constexpr size_t distNum  = 1;
    	ASSERT (fDict [labelNum] == "label");
    	ASSERT (fDict [distNum] == "dist");
    	
    	Steiner* parentNode = nullptr;
    	if (parent != no_id)
    	{
    	  parentNode = id2steiner [parent];
    	  if (! parentNode)
    	    throw runtime_error ("Parent id " + toString (parent) + " is not found");
    	}
    	
    	Real len = NaN;
    	if (parentNode && ! features [distNum]. empty ())
    	{
    	  len = max (0.0, str2real (features [distNum]));
    	  ASSERT (len >= 0.0);
    	}    	
    	
   	  const string name (features [labelNum]);
    	if (name. empty ())
    	  id2steiner [id] = new Steiner (tree, parentNode, len);
    	else
    	  tree. name2leaf [name] = new Leaf (tree, parentNode, len, name);
    }
};




struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Convert a tree in textual ASN.1 format to internal format")
	{
	  // Input
	  addPositional ("input_tree", "Tree in textual ASN.1 format (BioTreeContainer)");
	}



	void body () const final
  {
		const string input_tree = getArg ("input_tree");    


    DistTree tree;
    {
      TreeAsn asn (input_tree, tree);
      asn. asnRead ();
      const string title ("BioTreeContainer");
      if (asn. title != title)
        throw runtime_error ("Title " + strQuote (title) + " is expected");
    }
    tree. qc ();         
    tree. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


