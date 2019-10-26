// char2attr.cpp

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
*   Add a Boolean attribute for the object names containing a specified substring
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
    : Application ("Add a Boolean attribute for the object names containing a specified substring")
  	{
  	  version = VERSION;
  	  addFlag("equal_multiplicity", "the objects with the substring should have teh same weight as the other objects");
  	  addPositional ("in_file", "Input " + dmSuff + "-file");
  	  addPositional ("substr", "Object name substring");	
  	  addPositional ("attr_name", "Name of the new Boolean attribute");
  	  addPositional ("out_file", "Output " + dmSuff + "-file");
  	}
	
	
	
	void body () const final
	{
		const bool equal_multiplicity = getFlag ("equal_multiplicity");
		const string in_file          = getArg  ("in_file");
		const string substr           = getArg ("substr");
		const string attr_name        = getArg ("attr_name");
		const string out_file         = getArg ("out_file");
		ASSERT (! substr. empty ());
		ASSERT (! attr_name. empty ());


    Dataset ds (in_file);
    ExtBoolAttr1* attr = new ExtBoolAttr1 (attr_name, ds);   
    FOR (size_t, i, ds. objs. size ())
      (*attr) [i] = (ebool) contains (ds. objs. at (i) -> name, substr);
      
    if (equal_multiplicity)
    {
      Real classMult = 0;
      Real otherMult = 0;
      FOR (size_t, i, ds. objs. size ())
        if ((*attr) [i])
          classMult += ds. objs. at (i) -> mult;
        else
          otherMult += ds. objs. at (i) -> mult;
      ASSERT (classMult > 0);
      const Real ratio = otherMult / classMult;

      FOR (size_t, i, ds. objs. size ())
        if ((*attr) [i])
          const_cast <Obj*> (ds. objs. at (i)) -> mult *= ratio;
    }
      
    {
      OFStream of ("", out_file, dmExt);
      ds. saveText (of);
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


