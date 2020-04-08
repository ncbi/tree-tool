// feature2gain_loss.cpp

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
*   Print gains and losses of one feature in a tree
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "featureTree.hpp"
using namespace FeatureTree_sp;
#include "version.inc"



namespace 
{



struct ThisApplication : Application
{
	ThisApplication ()
  	: Application ("Print gains and losses of one feature in a tree")
  	{
  	  version = VERSION;
  	  addPositional ("input_tree", "Input file with the tree");
  	  addPositional ("feature", "Input file with genome list");
  	  addFlag ("prefer_gain", "Prefer gain over loss in maximum parsimony method");
  	}



	void body () const final
  {
		const string input_tree   = getArg ("input_tree");
		const string featureFName = getArg ("feature");
		const bool   prefer_gain  = getFlag ("prefer_gain");

		
    const FeatureTree tree (input_tree, featureFName, prefer_gain);
    tree. qc ();    

    ASSERT (tree. features. size () == 1);
    const Feature& feature = tree. features [0];
    for (const Phyl* phyl : feature. gains)
      cout << "gain\t" << phyl->getLcaName () << endl;
    if (feature. rootGain)
      cout << "gain\t" << tree. root->getLcaName () << endl;
    for (const Phyl* phyl : feature. losses)
      cout << "loss\t" << phyl->getLcaName () << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}


