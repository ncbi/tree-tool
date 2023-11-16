// dm_merge.cpp

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
*   Merge datasets
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Merge 2 " + dmSuff + "-files and print the resulting " + dmSuff + "-file. Files must have same objects in the same order, and different attributes")
  {
    version = VERSION;
	  addPositional ("file1", "First "  + dmSuff + "-file");
	  addPositional ("file2", "Second " + dmSuff + "-file");
	}



	void body () const final
	{
		const string fName1 = getArg ("file1");
		const string fName2 = getArg ("file2");
		
		
    Dataset ds1 (fName1);
    const Dataset ds2 (fName2);
    if (ds1. objs. size () != ds2. objs. size ())
    	throw runtime_error ("Datasets have different number of objects");
    FFOR (size_t, objNum, ds1. objs. size ())
      if (ds1. objs [objNum] -> name != ds2. objs [objNum] -> name)
	    	throw runtime_error ("Dataset 1 has object " + ds1. objs [objNum] -> name + ", but dataset 2 has object " + ds2. objs [objNum] -> name);

    for (const Attr* attr : ds2. attrs)
    	attr->copyToDataset (ds1);
    
    ds1. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



