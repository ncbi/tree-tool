// dm2objs.cpp

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
*   Print the list of objects of a dataset
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "version.inc"



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print the object names")
    {
      version = VERSION;
      addPositional ("file", dmSuff + "-file");
      addFlag ("comments", "Include comments");
    }
	
	
	
	void body () const
	{
		const string inFName = getArg ("file");
		const bool comments  = getFlag ("comments");


    const Dataset ds (inFName);
    
    for (const Obj* obj : ds. objs)
    {
    	cout << obj->name;
    	if (comments)
    	  cout << '\t' << obj->comment;
    	cout << endl;
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


