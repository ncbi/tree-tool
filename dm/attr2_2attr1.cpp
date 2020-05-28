// attr2_2attr1.cpp

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
*   Convert a two-way attribute into a one-way attribute
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
    : Application ("Print a " + dmSuff + "-file where all two-way attributes are converted into one-way attributes.\n\
Objects will be named: <objName1>-<objName2>")
    {
      version = VERSION;
  	  addPositional ("file", dmSuff + "-file");
  	  addFlag ("symmetrize", "Make all two-way attributes symmetric by averaging");
  	  addFlag ("zero_diagonal", "Skip diagonal values because they are 0");
  	}



	void body () const final
	{
		const string fName       = getArg ("file");
		const bool symmetrize    = getFlag ("symmetrize");
		const bool zero_diagonal = getFlag ("zero_diagonal");
		
		
    const Dataset ds (fName);
    
    Dataset ds1;
    FFOR (size_t, i, ds. objs. size ())
      FOR (size_t, j, symmetrize ? i + 1 : ds. objs. size ())
        if (! zero_diagonal || i != j)
        {
          const size_t objNum = ds1. appendObj (ds. pair2name (i, j, symmetrize));
          var_cast (ds1. objs [objNum]) -> mult =   ds. objs [i] -> mult
                                                  * ds. objs [j] -> mult;
        }
    ds1. setName2objNum ();
    
    for (const Attr* attr : ds. attrs)
      if (const Attr2* attr2 = attr->asAttr2 ())
      {
        Attr1* attr1 = attr2->createAttr1 (ds1);
        if (symmetrize)
          var_cast (attr2) -> symmetrize ();
        FFOR (size_t, i, ds. objs. size ())
          FFOR (size_t, j, symmetrize ? i + 1 : ds. objs. size ())
            if (! zero_diagonal || i != j)
              {
                const size_t objNum = ds1. getName2objNum (ds. pair2name (i, j, symmetrize));
                ASSERT (objNum != NO_INDEX);
                attr1->str2value (objNum, attr2->value2str (i, j));
              }
      }
      
    ds1. qc ();
    
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



