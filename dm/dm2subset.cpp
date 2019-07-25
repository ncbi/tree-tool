// dm2subset.cpp

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
*   Make a subset of objects
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a subset of a " + dmSuff + "-file")
  {
	  addFlag ("exclude", "Exclude the list of objects, otherwise include");
	  addPositional ("file", dmSuff + "-file");
	  addPositional ("objNameFName", "File with object names defining the subset");
	}



	void body () const
	{
		const bool   exclude      = getFlag ("exclude");
		const string fName        = getArg ("file");
		const string objNameFName = getArg ("objNameFName");
		
		
    const Dataset ds (fName);
    
    Set<string> objNames;
    {
      LineInput li (objNameFName);
      while (li. nextLine ())
        objNames << li. line;
    }
    
    Sample sm (ds);
    FOR (size_t, row, ds. objs. size ())
      if (objNames. contains (ds. objs [row] -> name) == exclude)
        sm. mult [row] = 0;
    sm. finish ();    

    sm. save (VectorPtr<Attr> (ds. attrs), cout);  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



