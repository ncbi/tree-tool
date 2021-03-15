// hmmNAME2ACC.cpp

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
*   Add ACC=NAME in an HMM library
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Add ACC=NAME in an HMM library")
  	{
      version = VERSION;
  	  addPositional ("in", "Input NCBI HMM library");
  	  addPositional ("out", "Output NCBI HMM library");
  	}



  void body () const final
	{
		const string in  = getArg ("in");
		const string out = getArg ("out");


    OFStream ofs (out);
    LineInput li (in);
    Progress prog (0, 1000);  // PAR
    while (li. nextLine ())
    {
      string s (li. line);
   	  replace (s, '\t', ' ');
      if (isLeft (s, "ACC "))
        continue; 
      if (isLeft (s, "NAME "))
      {
        prog ();
        findSplit (s);
        trim (s);
        ofs << "ACC\t" << s << endl;
      }
      ofs << li. line << endl;
    }
  }
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



