// randomDistTree.cpp

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
*   Make a random distance tree
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


struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Print a random distance tree")
	{
	  version = VERSION;
	  // Input
	  addPositional ("branch_prob", "Probability to expand a branch");
	  addPositional ("leaf_num_max", "Max. number of leaves");
	}



	void body () const final
  {
		const Real branch_prob    = str2<Prob> (getArg ("branch_prob"));
		const size_t leaf_num_max = str2<size_t> (getArg ("leaf_num_max"));
		ASSERT (isProb (branch_prob));
		ASSERT (branch_prob < 1.0);
		ASSERT (branch_prob > 0.0);
		ASSERT ((bool) leaf_num_max);
    

    DistTree tree (branch_prob, leaf_num_max);
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


