// dm2space.cpp

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
*   Make a subset of attributes
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a space of a " + dmSuff + "-file")
  {
    version = VERSION;
	  addPositional ("file", dmSuff + "-file");
	  addPositional ("attrNameFName", "File with attribute names defining the space");
	  addFlag ("exclude", "Exclude the list of attributes, otherwise include");
	}



	void body () const final
	{
		const string fName         = getArg ("file");
		const string attrNameFName = getArg ("attrNameFName");
		const bool   exclude       = getFlag ("exclude");
		
				
    Set<string> attrNames;  // This is faster than Dataset()
    {
      LineInput li (attrNameFName);
      while (li. nextLine ())
      {
      	trim (li. line);
      	if (! li. line. empty ())
          attrNames << li. line;
      }
    }
    
    const Dataset ds (fName);
    
    VectorPtr<Attr> attrs;  attrs. reserve (ds. attrs. size ());
    for (const string& name : attrNames)
    	if (const Attr* attr = ds. name2attr (name))
    		attrs << attr;
    	else
    	  throw runtime_error ("Attribute " + name + " is not in the dataset");
    	  
    if (exclude)
    {
    	VectorPtr<Attr> attrs1 (ds. attrs);
    	attrs. sort ();
    	attrs1. setMinus (attrs);
    	attrs = attrs1;
    }
    
    const Sample sm (ds);
    sm. save (nullptr, attrs, cout);  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



