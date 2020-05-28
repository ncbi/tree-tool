// normalCdf.cpp

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
*   C.d.f. of a normal distribution
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
    : Application ("Normal distribution CDF")
  	{
  	  version = VERSION;
  	  addPositional ("loc",   "Location parameter");	
  	  addPositional ("scale", "Scale parameter");	  
  	  addPositional ("x",     "Value");
  	}



	void body () const final
	{
		const Real loc   = str2real (getArg ("loc"));
		const Real scale = str2real (getArg ("scale"));
		const Real x     = str2real (getArg ("x"));

    Normal norm;
    norm. setParam (loc, scale);
    cout << norm. cdf (x) << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


