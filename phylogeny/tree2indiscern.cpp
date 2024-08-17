// tree2indiscern.cpp

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
*   Print indiscernibility clusters of a distance tree
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


const string distName = "dist"; 



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Print indiscernibility clusters of a distance tree: <object> <cluster leading object>")
	{
	  version = VERSION;
	  addPositional ("input_tree", "File with the tree and arc lengths");
	  addFlag ("discern", "Print discernible objects");
	}



	void body () const final
  {
		const string input_tree  = getArg ("input_tree");
		const bool   discernible = getFlag ("discern");
				
    
    const DistTree tree (input_tree); 
    tree. qc (); 
    
    Set<string> repr;
    for (const auto& it : tree. name2leaf)
    {
      const Leaf* leaf = it. second;
      if (discernible)
      {
        if (leaf->discernible)
          cout << leaf->name << endl;
        else
        {
          string s (static_cast <const DTNode*> (leaf->getParent ()) -> getFirstDecendant () -> getName ());
          if (! repr. contains (s))
          {
            cout << s << endl;
            repr << std::move (s);
          }
        }
      }
      else
        if (! leaf->discernible)
        {
          const string s1 (leaf->name);
          const string s2 (static_cast <const DTNode*> (leaf->getParent ()) -> getFirstDecendant () -> getName ());
          cout << s1 << '\t' << s2 << endl;
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


