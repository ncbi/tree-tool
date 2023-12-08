// unCgi.cpp

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
*   Extract URL parameters
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"

#include "cgi.hpp"
using namespace Cgi_sp;

#include "../common.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Save HTTP parameters as files <prefix>_<parameter>")
	  // Why not <prefix>.<parameter> ??
  	{
      version = VERSION;
  	  addPositional ("prefix", "string");
  	}



	void body () const final
  {
		const string prefix = getArg ("prefix");
		if (verbose ())
		  cout << prefix << endl;
  

    const KeyValue keyValues (unCgi ());
    for (const auto& it : keyValues)
    {
      QC_ASSERT (! it. first. empty ());
      OFStream f (prefix + "_" + it. first);
      f << it. second << endl;
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


