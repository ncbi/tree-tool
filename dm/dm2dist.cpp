// pca.cpp

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
*   Create a 2-way distance attribute
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace
{



struct ThisApplication : Application
{	
  ThisApplication ()
  : Application ("Create a 2-way distance attribute")
	{
	  version = VERSION;
	  
	  addPositional ("file", dmSuff + "-file without the extension");	  
	  addPositional ("dist", "Name of created 2-way attribute of squared distance in L_2");
	  addFlag ("standardize", "Standardize the attributes");
	}
	
	
	
	void body () const final
	{
		const string inFName      = getArg("file");
		const string distAttrName = getArg("dist");
		const bool   stndP        = getFlag ("standardize");
		
    
    Dataset ds (inFName);
    ds. qc ();
    
    if (ds. attrs. size () <= 0)
      throw runtime_error ("No attributes");
    
    const Sample sm_orig (ds);
    Sample sm (ds);
    if (sm. mult_sum <= 0.0)
      throw runtime_error ("Too small data size");

    Space1<Attr1> spRaw (ds, true);
    if (spRaw. empty ())
      throw runtime_error ("No attributes");
    spRaw. qc ();
    
    Space1<NumAttr1> sp;
    if (stndP)
      sp = std::move (spRaw. standardize (sm, ds));  
    else
      sp = std::move (spRaw. toNumAttr1 (ds));
        
    sp. qc ();
      
    getDist2 (sp, distAttrName, ds);
    ds. qc ();

    ds. saveText (cout);    
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


